use crate::config::{ConsensusMode, PoaConfig};
use crate::error::PoaError;
use crate::types::{Consensus, GraphStats};
use std::collections::HashMap;

// ─── Sentinel ────────────────────────────────────────────────────────────────

/// "Not yet filled" sentinel. Safe to add gap_extend once without i32 overflow.
const UNSET: i32 = i32::MIN / 2;

#[inline]
fn safe_add(a: i32, b: i32) -> i32 {
    if a == UNSET {
        UNSET
    } else {
        a.saturating_add(b)
    }
}

/// Sentinel topo index meaning "came from the virtual start node".
const VIRTUAL: usize = usize::MAX;

// ─── Band helpers ─────────────────────────────────────────────────────────────

/// Returned by `compute_effective_band` to mean "no band restriction".
const UNBANDED: usize = usize::MAX;

fn compute_effective_band(cfg: &PoaConfig, read_len: usize) -> usize {
    if cfg.adaptive_band {
        let w = cfg.adaptive_band_b + (cfg.adaptive_band_f * read_len as f32).ceil() as usize;
        let w = if cfg.band_width > 0 {
            w.max(cfg.band_width)
        } else {
            w
        };
        w.max(1)
    } else if cfg.band_width > 0 {
        cfg.band_width
    } else {
        UNBANDED
    }
}

/// Lowest j>=1 reachable from topo position `t` within band `w`.
/// The j=0 column is handled separately, gated via `in_band_at`.
#[inline]
fn j_lo_at(t: usize, w: usize) -> usize {
    if w == UNBANDED {
        1
    } else {
        t.saturating_sub(w).max(1)
    }
}

/// Highest j reachable from topo position `t` within band `w`.
#[inline]
fn j_hi_at(t: usize, w: usize, l: usize) -> usize {
    if w == UNBANDED { l } else { (t + w).min(l) }
}

/// True if cell (t, j) is within the diagonal band.
/// j=0 is in-band only when `t <= w` (i.e. `t.saturating_sub(w) == 0`).
#[inline]
fn in_band_at(t: usize, j: usize, w: usize, l: usize) -> bool {
    if w == UNBANDED {
        return true;
    }
    j >= t.saturating_sub(w) && j <= (t + w).min(l)
}

/// Cells per row in band-indexed storage.
/// Unbanded: l+1 (direct j indexing). Banded: 2*w+2 (index 0 = j=0, index k>=1 = j=j_lo+k-1).
#[inline]
fn row_stride_for(w: usize, l: usize) -> usize {
    if w == UNBANDED { l + 1 } else { 2 * w + 2 }
}

/// Local column index within a band-indexed row.
/// Unbanded: returns `j` directly. Banded: j=0 -> 0, j>=1 -> 1 + (j - j_lo_at(t, w)).
#[inline]
fn local_j(t: usize, j: usize, w: usize) -> usize {
    if w == UNBANDED {
        j
    } else if j == 0 {
        0
    } else {
        j - j_lo_at(t, w) + 1
    }
}

/// DP table backed by a flat Vec with band-indexed row layout.
struct DpTable {
    data: Vec<Cell>,
    rs: usize,
    w: usize,
}

impl DpTable {
    fn new(n: usize, w: usize, l: usize) -> Self {
        let rs = row_stride_for(w, l);
        DpTable {
            data: vec![Cell::unset(); n * rs],
            rs,
            w,
        }
    }

    #[inline]
    fn get(&self, t: usize, j: usize) -> Cell {
        self.data[t * self.rs + local_j(t, j, self.w)]
    }

    #[inline]
    fn set(&mut self, t: usize, j: usize, cell: Cell) {
        let idx = t * self.rs + local_j(t, j, self.w);
        self.data[idx] = cell;
    }
}

// ─── Graph types ─────────────────────────────────────────────────────────────

struct Node {
    base: u8,
    /// (successor_node_index, edge_weight)
    out_edges: Vec<(usize, i32)>,
    /// predecessor node indices
    in_edges: Vec<usize>,
    /// reads that produced a Match op at this node (not Delete ops)
    coverage: u32,
    /// reads that produced a Delete op at this node (traversed without consuming a base)
    delete_count: u32,
}

pub struct PoaGraph {
    nodes: Vec<Node>,
    config: PoaConfig,
    /// number of reads integrated (seed + subsequent)
    n_reads: usize,
    /// original reads stored for per-allele sub-graph reconstruction
    reads: Vec<Vec<u8>>,
    /// (from_node, to_node) → sorted list of read indices that traversed that edge
    edge_reads: HashMap<(usize, usize), Vec<u32>>,
    /// number of times the long-unbanded warning was emitted
    warnings: usize,
}

// ─── DP cell types ───────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq, Eq)]
enum State {
    M,
    I,
    D,
}

#[derive(Clone, Copy)]
struct Cell {
    score: i32,
    /// topo rank of predecessor node, or VIRTUAL
    pred_t: usize,
}

impl Cell {
    fn unset() -> Self {
        Cell {
            score: UNSET,
            pred_t: 0,
        }
    }
}

/// Alignment operation produced by traceback, in forward read order.
#[derive(Clone, Copy, Debug)]
pub enum AlignOp {
    /// Base matched or mismatched against this graph node index.
    Match(usize),
    /// Inserted base from the read, no corresponding graph node.
    Insert(u8),
    /// Graph node skipped by the read (gap in read).
    Delete(usize),
}

// ─── Graph construction helpers ──────────────────────────────────────────────

fn push_node(nodes: &mut Vec<Node>, base: u8) -> usize {
    let idx = nodes.len();
    nodes.push(Node {
        base,
        out_edges: Vec::new(),
        in_edges: Vec::new(),
        coverage: 0,
        delete_count: 0,
    });
    idx
}

fn add_edge(nodes: &mut [Node], from: usize, to: usize) {
    nodes[from].out_edges.push((to, 1));
    nodes[to].in_edges.push(from);
}

fn increment_or_add_edge(nodes: &mut [Node], from: usize, to: usize) {
    for (succ, w) in nodes[from].out_edges.iter_mut() {
        if *succ == to {
            *w += 1;
            return;
        }
    }
    nodes[from].out_edges.push((to, 1));
    nodes[to].in_edges.push(from);
}

// ─── Topological sort ────────────────────────────────────────────────────────

fn topological_order(nodes: &[Node]) -> (Vec<usize>, Vec<usize>) {
    let n = nodes.len();
    let mut in_deg: Vec<usize> = nodes.iter().map(|nd| nd.in_edges.len()).collect();
    let mut queue: std::collections::VecDeque<usize> = (0..n).filter(|&i| in_deg[i] == 0).collect();
    let mut topo: Vec<usize> = Vec::with_capacity(n);

    while let Some(u) = queue.pop_front() {
        topo.push(u);
        for &(v, _) in &nodes[u].out_edges {
            in_deg[v] -= 1;
            if in_deg[v] == 0 {
                queue.push_back(v);
            }
        }
    }

    let mut rank_of = vec![0usize; n];
    for (t, &node_idx) in topo.iter().enumerate() {
        rank_of[node_idx] = t;
    }
    (topo, rank_of)
}

// ─── DP alignment ────────────────────────────────────────────────────────────

/// Align `query` into the graph using band width `w` (`UNBANDED` = full NW).
/// Returns alignment ops in forward order, or `Err(BandTooNarrow)` when the
/// traceback path reaches the band edge (alignment result would be incorrect).
fn align(
    nodes: &[Node],
    topo: &[usize],
    rank_of: &[usize],
    query: &[u8],
    cfg: &PoaConfig,
    w: usize,
) -> Result<Vec<AlignOp>, PoaError> {
    let n = topo.len();
    let l = query.len();

    let go = cfg.gap_open; // negative
    let ge = cfg.gap_extend; // negative

    let mut m = DpTable::new(n, w, l);
    let mut ins = DpTable::new(n, w, l);
    let mut del = DpTable::new(n, w, l);

    for (t, &node_idx) in topo.iter().enumerate() {
        let node_base = nodes[node_idx].base;
        let is_source = nodes[node_idx].in_edges.is_empty();

        // ── j = 0: del only; in-band when t <= w ──────────────────────────────
        if in_band_at(t, 0, w, l) {
            let (mut best, mut best_pred) = (UNSET, 0usize);
            if is_source {
                let val = go + ge;
                if val > best {
                    best = val;
                    best_pred = VIRTUAL;
                }
            }
            for &p in &nodes[node_idx].in_edges {
                let p_t = rank_of[p];
                let vm = safe_add(m.get(p_t, 0).score, go + ge);
                if vm > best {
                    best = vm;
                    best_pred = p_t;
                }
                let vd = safe_add(del.get(p_t, 0).score, ge);
                if vd > best {
                    best = vd;
                    best_pred = p_t;
                }
            }
            if best != UNSET {
                del.set(
                    t,
                    0,
                    Cell {
                        score: best,
                        pred_t: best_pred,
                    },
                );
            }
        }

        // ── j = j_lo..=j_hi ───────────────────────────────────────────────────
        let j_lo = j_lo_at(t, w);
        let j_hi = j_hi_at(t, w, l);
        for j in j_lo..=j_hi {
            let q_base = query[j - 1];
            let score_fn = if node_base == q_base {
                cfg.match_score
            } else {
                cfg.mismatch_score
            };

            // M[t][j] ─────────────────────────────────────────────────────────
            {
                let (mut best, mut best_pred) = (UNSET, VIRTUAL);
                if is_source && j == 1 {
                    let val = score_fn;
                    if val > best {
                        best = val;
                        best_pred = VIRTUAL;
                    }
                }
                if is_source && j > 1 {
                    let val = safe_add(go + (j as i32 - 1) * ge, score_fn);
                    if val > best {
                        best = val;
                        best_pred = VIRTUAL;
                    }
                }
                for &p in &nodes[node_idx].in_edges {
                    let p_t = rank_of[p];
                    if in_band_at(p_t, j - 1, w, l) {
                        let vm = safe_add(m.get(p_t, j - 1).score, score_fn);
                        if vm > best {
                            best = vm;
                            best_pred = p_t;
                        }
                        let vi = safe_add(ins.get(p_t, j - 1).score, score_fn);
                        if vi > best {
                            best = vi;
                            best_pred = p_t;
                        }
                        let vd = safe_add(del.get(p_t, j - 1).score, score_fn);
                        if vd > best {
                            best = vd;
                            best_pred = p_t;
                        }
                    }
                }
                if best != UNSET {
                    m.set(
                        t,
                        j,
                        Cell {
                            score: best,
                            pred_t: best_pred,
                        },
                    );
                }
            }

            // I[t][j] ─────────────────────────────────────────────────────────
            {
                let (mut best, mut best_pred) = (UNSET, VIRTUAL);
                if is_source && j == 1 {
                    let val = go + ge;
                    if val > best {
                        best = val;
                        best_pred = VIRTUAL;
                    }
                }
                if in_band_at(t, j - 1, w, l) {
                    let vm = safe_add(m.get(t, j - 1).score, go + ge);
                    if vm > best {
                        best = vm;
                        best_pred = t;
                    }
                    let vi = safe_add(ins.get(t, j - 1).score, ge);
                    if vi > best {
                        best = vi;
                        best_pred = t;
                    }
                }
                if best != UNSET {
                    ins.set(
                        t,
                        j,
                        Cell {
                            score: best,
                            pred_t: best_pred,
                        },
                    );
                }
            }

            // D[t][j] ─────────────────────────────────────────────────────────
            // No is_source init here: source-deletion initialization belongs only
            // in the j=0 column. Setting D[source][j≥1]=const here would let the
            // traceback reach D[source] at any j>0 and break early via pred=VIRTUAL,
            // leaving read characters unconsumed and producing wrong alignment ops.
            {
                let (mut best, mut best_pred) = (UNSET, 0usize);
                for &p in &nodes[node_idx].in_edges {
                    let p_t = rank_of[p];
                    if in_band_at(p_t, j, w, l) {
                        let vm = safe_add(m.get(p_t, j).score, go + ge);
                        if vm > best {
                            best = vm;
                            best_pred = p_t;
                        }
                        let vi = safe_add(ins.get(p_t, j).score, go + ge);
                        if vi > best {
                            best = vi;
                            best_pred = p_t;
                        }
                        let vd = safe_add(del.get(p_t, j).score, ge);
                        if vd > best {
                            best = vd;
                            best_pred = p_t;
                        }
                    }
                }
                if best != UNSET {
                    del.set(
                        t,
                        j,
                        Cell {
                            score: best,
                            pred_t: best_pred,
                        },
                    );
                }
            }
        }
    }

    // ── Find best terminal cell at column l ───────────────────────────────────
    // Both Global and SemiGlobal scan all topo positions at column l.
    // (SemiGlobal free-end-gaps are achieved by not penalising leading/trailing
    // graph nodes that are skipped; the DP initialisation handles free start and
    // this terminal scan handles free end.)
    // When every topo position has j=l outside its band, the alignment cannot
    // complete — the band is definitively too narrow.
    let terminal_best = (0..n)
        .filter(|&t| in_band_at(t, l, w, l))
        .flat_map(|t| {
            [
                (t, State::M, m.get(t, l).score),
                (t, State::I, ins.get(t, l).score),
                (t, State::D, del.get(t, l).score),
            ]
        })
        .filter(|&(_, _, sc)| sc != UNSET)
        .max_by_key(|&(_, _, sc)| sc)
        .map(|(t, s, _)| (t, s));

    let (best_t, best_state) = match terminal_best {
        Some(result) => result,
        None if w != UNBANDED && l > 0 => {
            return Err(PoaError::BandTooNarrow {
                band_width: w,
                read_len: l,
            });
        }
        None => (0, State::M),
    };

    // ── Traceback ─────────────────────────────────────────────────────────────
    let mut ops: Vec<AlignOp> = Vec::with_capacity(l + n / 4);
    let mut t = best_t;
    let mut j = l;
    let mut cur_state = best_state;

    loop {
        let cell = match cur_state {
            State::M => m.get(t, j),
            State::I => ins.get(t, j),
            State::D => del.get(t, j),
        };

        if cell.score == UNSET {
            break;
        }

        match cur_state {
            State::M => {
                ops.push(AlignOp::Match(topo[t]));
                if cell.pred_t == VIRTUAL {
                    for k in (1..j).rev() {
                        ops.push(AlignOp::Insert(query[k - 1]));
                    }
                    break;
                }
                t = cell.pred_t;
                j -= 1;
                cur_state = best_prev_state(&m, &ins, &del, t, j);
            }
            State::I => {
                ops.push(AlignOp::Insert(query[j - 1]));
                j -= 1;
                if cell.pred_t == VIRTUAL {
                    for k in (1..j).rev() {
                        ops.push(AlignOp::Insert(query[k - 1]));
                    }
                    break;
                }
                cur_state = if m.get(t, j).score >= ins.get(t, j).score {
                    State::M
                } else {
                    State::I
                };
            }
            State::D => {
                ops.push(AlignOp::Delete(topo[t]));
                if cell.pred_t == VIRTUAL {
                    break;
                }
                t = cell.pred_t;
                cur_state = best_prev_state(&m, &ins, &del, t, j);
            }
        }

        if j == 0 && cur_state != State::D {
            break;
        }
    }

    ops.reverse();

    Ok(ops)
}

#[inline]
fn best_prev_state(m: &DpTable, ins: &DpTable, del: &DpTable, t: usize, j: usize) -> State {
    let sm = m.get(t, j).score;
    let si = ins.get(t, j).score;
    let sd = del.get(t, j).score;
    if sm != UNSET && sm >= si && sm >= sd {
        State::M
    } else if si != UNSET && si >= sd {
        State::I
    } else {
        State::D
    }
}

// ─── Graph update ─────────────────────────────────────────────────────────────

fn add_to_graph(
    nodes: &mut Vec<Node>,
    edge_reads: &mut HashMap<(usize, usize), Vec<u32>>,
    query: &[u8],
    ops: &[AlignOp],
    read_idx: u32,
) {
    let mut prev: Option<usize> = None;
    let mut q_idx: usize = 0;

    for &op in ops {
        match op {
            AlignOp::Match(node_idx) => {
                let q_base = query[q_idx];
                q_idx += 1;
                let cur = if nodes[node_idx].base == q_base {
                    nodes[node_idx].coverage += 1;
                    if let Some(p) = prev {
                        increment_or_add_edge(nodes, p, node_idx);
                        edge_reads.entry((p, node_idx)).or_default().push(read_idx);
                    }
                    node_idx
                } else {
                    let new_idx = push_node(nodes, q_base);
                    nodes[new_idx].coverage = 1;
                    if let Some(p) = prev {
                        nodes[p].out_edges.push((new_idx, 1));
                        nodes[new_idx].in_edges.push(p);
                        edge_reads.entry((p, new_idx)).or_default().push(read_idx);
                    }
                    new_idx
                };
                prev = Some(cur);
            }
            AlignOp::Insert(q_base) => {
                q_idx += 1;
                let new_idx = push_node(nodes, q_base);
                nodes[new_idx].coverage = 1;
                if let Some(p) = prev {
                    nodes[p].out_edges.push((new_idx, 1));
                    nodes[new_idx].in_edges.push(p);
                    edge_reads.entry((p, new_idx)).or_default().push(read_idx);
                }
                prev = Some(new_idx);
            }
            AlignOp::Delete(node_idx) => {
                // Traverse without consuming a query base; counts as a gap vote for MF.
                nodes[node_idx].delete_count += 1;
                if let Some(p) = prev {
                    increment_or_add_edge(nodes, p, node_idx);
                    edge_reads.entry((p, node_idx)).or_default().push(read_idx);
                }
                prev = Some(node_idx);
            }
        }
    }
}

// ─── Heaviest path ────────────────────────────────────────────────────────────

fn heaviest_path(nodes: &[Node], topo: &[usize], rank_of: &[usize]) -> Vec<(usize, u8)> {
    let n = topo.len();
    let mut cum: Vec<(i64, Option<usize>)> = vec![(0, None); n];

    for t in 0..n {
        let node_idx = topo[t];
        let curr = cum[t].0;
        for &(succ_idx, weight) in &nodes[node_idx].out_edges {
            let succ_t = rank_of[succ_idx];
            let candidate = curr + (weight - 1) as i64;
            if candidate > cum[succ_t].0 {
                cum[succ_t] = (candidate, Some(t));
            }
        }
    }

    let max_cum = (0..n).map(|t| cum[t].0).max().unwrap_or(0);
    // Prefer shortest path among equal-weight termini.
    let best_t = (0..n).find(|&t| cum[t].0 == max_cum).unwrap_or(0);

    let mut path: Vec<(usize, u8)> = Vec::new();
    let mut t = best_t;
    loop {
        path.push((topo[t], nodes[topo[t]].base));
        match cum[t].1 {
            None => break,
            Some(pred_t) => t = pred_t,
        }
    }
    path.reverse();
    path
}

// ─── Majority-frequency consensus ────────────────────────────────────────────

/// Column-wise plurality vote: include a node when its base votes outnumber its
/// gap votes AND enough reads voted at that column to meet `min_cov`.
///
/// `coverage` = reads that matched the node. `delete_count` = reads that deleted
/// it (traversed without consuming a base). Reads that took a different branch
/// through a bubble don't vote at this column at all, so the denominator is
/// `coverage + delete_count`, not `n_reads`.
fn majority_frequency(nodes: &[Node], topo: &[usize], min_cov: u32) -> Vec<(usize, u8)> {
    topo.iter()
        .copied()
        .filter(|&idx| {
            let cov = nodes[idx].coverage;
            let del = nodes[idx].delete_count;
            let total = cov + del;
            total >= min_cov && cov * 2 >= total
        })
        .map(|idx| (idx, nodes[idx].base))
        .collect()
}

// ─── Graph statistics ─────────────────────────────────────────────────────────

fn compute_stats(nodes: &[Node], min_allele_freq: f64, n_reads: usize) -> GraphStats {
    let node_count = nodes.len();
    let edge_count: usize = nodes.iter().map(|nd| nd.out_edges.len()).sum();

    // Coverage mean and variance
    let coverages: Vec<f64> = nodes.iter().map(|nd| nd.coverage as f64).collect();
    let coverage_mean = if node_count == 0 {
        0.0
    } else {
        coverages.iter().sum::<f64>() / node_count as f64
    };
    let coverage_variance = if node_count == 0 {
        0.0
    } else {
        coverages
            .iter()
            .map(|&c| (c - coverage_mean).powi(2))
            .sum::<f64>()
            / node_count as f64
    };

    // Single-support fraction
    let single_support = nodes.iter().filter(|nd| nd.coverage == 1).count();
    let single_support_fraction = if node_count == 0 {
        0.0
    } else {
        single_support as f64 / node_count as f64
    };

    // Edge weight Gini coefficient
    let mut weights: Vec<f64> = nodes
        .iter()
        .flat_map(|nd| nd.out_edges.iter().map(|&(_, w)| w as f64))
        .collect();
    weights.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let edge_weight_gini = if weights.len() < 2 {
        0.0
    } else {
        let n = weights.len() as f64;
        let sum: f64 = weights.iter().sum();
        if sum == 0.0 {
            0.0
        } else {
            let numerator: f64 = weights
                .iter()
                .enumerate()
                .map(|(i, &w)| (2.0 * (i as f64 + 1.0) - n - 1.0) * w)
                .sum::<f64>();
            numerator / (n * sum)
        }
    };

    // Bubble count and max_bubble_depth.
    // A bubble entry is a node with 2+ out-edges each above the min_allele_freq threshold.
    // max_bubble_depth = the highest minority-arm weight across all bubbles, i.e. the
    // strongest signal for a second allele.
    let threshold = (n_reads as f64 * min_allele_freq).ceil() as i32;
    let mut bubble_count = 0usize;
    let mut max_bubble_depth = 0usize;
    for nd in nodes {
        let mut qualifying: Vec<i32> = nd
            .out_edges
            .iter()
            .filter(|&&(_, w)| w >= threshold)
            .map(|&(_, w)| w)
            .collect();
        if qualifying.len() >= 2 {
            bubble_count += 1;
            qualifying.sort_unstable_by(|a, b| b.cmp(a)); // descending
            max_bubble_depth = max_bubble_depth.max(qualifying[1] as usize);
        }
    }

    // Mean per-column Shannon entropy (bits).
    // At each node, the "column" has two outcomes: base (coverage) or gap (delete_count).
    // Entropy is 0 when all reads agree and 1 bit when perfectly split.
    // Nodes that were never voted on (coverage=0, delete_count=0) are skipped.
    let mean_column_entropy = {
        let mut sum = 0.0f64;
        let mut count = 0usize;
        for nd in nodes {
            let cov = nd.coverage as f64;
            let del = nd.delete_count as f64;
            let total = cov + del;
            if total > 0.0 {
                let h = binary_entropy(cov / total);
                sum += h;
                count += 1;
            }
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    };

    GraphStats {
        node_count,
        edge_count,
        bubble_count,
        max_bubble_depth,
        coverage_mean,
        coverage_variance,
        edge_weight_gini,
        single_support_fraction,
        mean_column_entropy,
    }
}

#[inline]
fn binary_entropy(p: f64) -> f64 {
    if p <= 0.0 || p >= 1.0 {
        0.0
    } else {
        let q = 1.0 - p;
        -(p * p.log2() + q * q.log2())
    }
}

// ─── Multi-allele helpers ─────────────────────────────────────────────────────

/// Returns `(entry_node_idx, arm_start_node_idxs)` for every bubble in the
/// graph — nodes that have 2+ outgoing edges each above the min_allele_freq
/// threshold. The list is in topological order.
fn find_bubbles(
    nodes: &[Node],
    topo: &[usize],
    n_reads: usize,
    min_allele_freq: f64,
) -> Vec<(usize, Vec<usize>)> {
    let threshold = ((n_reads as f64 * min_allele_freq).ceil() as i32).max(1);
    topo.iter()
        .copied()
        .filter_map(|node_idx| {
            let arms: Vec<usize> = nodes[node_idx]
                .out_edges
                .iter()
                .filter(|&&(_, w)| w >= threshold)
                .map(|&(to, _)| to)
                .collect();
            if arms.len() >= 2 {
                Some((node_idx, arms))
            } else {
                None
            }
        })
        .collect()
}

/// Partition reads into allele groups by examining which arm of `entry_node`
/// each read traversed, using per-edge read membership recorded in `edge_reads`.
/// Reads not observed on any qualifying arm are absorbed into the largest group.
/// Returns at most 2 groups (the two highest-support arms).
fn partition_reads_by_bubble(
    edge_reads: &HashMap<(usize, usize), Vec<u32>>,
    entry_node: usize,
    arm_starts: &[usize],
    n_reads: usize,
) -> Vec<Vec<usize>> {
    // Keep the two arms with the most read support.
    let mut arm_order: Vec<usize> = (0..arm_starts.len()).collect();
    arm_order.sort_unstable_by(|&a, &b| {
        let wa = edge_reads
            .get(&(entry_node, arm_starts[a]))
            .map_or(0, |v| v.len());
        let wb = edge_reads
            .get(&(entry_node, arm_starts[b]))
            .map_or(0, |v| v.len());
        wb.cmp(&wa)
    });
    let n_arms = arm_order.len().min(2);

    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); n_arms];
    let mut assigned = vec![false; n_reads];

    for (slot, &arm_idx) in arm_order[..n_arms].iter().enumerate() {
        let arm_start = arm_starts[arm_idx];
        if let Some(reads) = edge_reads.get(&(entry_node, arm_start)) {
            for &r in reads {
                let r = r as usize;
                if r < n_reads && !assigned[r] {
                    groups[slot].push(r);
                    assigned[r] = true;
                }
            }
        }
    }

    // Reads not on any arm (e.g. Delete traversals that skipped the entry)
    // go to the largest group.
    let largest = groups
        .iter()
        .enumerate()
        .max_by_key(|(_, g)| g.len())
        .map(|(i, _)| i)
        .unwrap_or(0);
    for (r, &done) in assigned.iter().enumerate() {
        if !done {
            groups[largest].push(r);
        }
    }

    groups.retain(|g| !g.is_empty());
    groups
}

/// Index within `group` of the read whose length is closest to the median.
fn choose_seed(group: &[usize], reads: &[Vec<u8>]) -> usize {
    if group.len() == 1 {
        return 0;
    }
    let mut lens: Vec<usize> = group.iter().map(|&i| reads[i].len()).collect();
    lens.sort_unstable();
    let median = lens[lens.len() / 2];
    group
        .iter()
        .enumerate()
        .min_by_key(|&(_, &i)| reads[i].len().abs_diff(median))
        .map(|(slot, _)| slot)
        .unwrap_or(0)
}

// ─── PoaGraph public API ─────────────────────────────────────────────────────

impl PoaGraph {
    pub fn new(seed: &[u8], config: PoaConfig) -> Result<Self, PoaError> {
        if seed.is_empty() {
            return Err(PoaError::EmptyInput);
        }

        let mut nodes: Vec<Node> = Vec::with_capacity(seed.len());
        for &base in seed {
            let idx = push_node(&mut nodes, base);
            nodes[idx].coverage = 1;
        }
        let n = nodes.len();
        // Record the seed (read 0) as traversing each backbone edge.
        let mut edge_reads: HashMap<(usize, usize), Vec<u32>> = HashMap::new();
        for i in 0..n.saturating_sub(1) {
            add_edge(&mut nodes, i, i + 1);
            edge_reads.insert((i, i + 1), vec![0u32]);
        }

        Ok(PoaGraph {
            nodes,
            config,
            n_reads: 1,
            reads: vec![seed.to_vec()],
            edge_reads,
            warnings: 0,
        })
    }

    pub fn add_read(&mut self, read: &[u8]) -> Result<(), PoaError> {
        if read.is_empty() {
            return Err(PoaError::EmptyInput);
        }

        let w = compute_effective_band(&self.config, read.len());

        if self.config.warn_on_long_unbanded && w == UNBANDED && read.len() > 1000 {
            eprintln!(
                "poa-consensus: warning: read length {} bp with band_width=0 (unbanded). \
                 Memory usage is O(n*m). Consider setting band_width or adaptive_band=true. \
                 Suppress with warn_on_long_unbanded=false.",
                read.len()
            );
            self.warnings += 1;
        }

        let (topo, rank_of) = topological_order(&self.nodes);
        let ops = align(&self.nodes, &topo, &rank_of, read, &self.config, w)?;
        let read_idx = self.n_reads as u32;
        add_to_graph(&mut self.nodes, &mut self.edge_reads, read, &ops, read_idx);
        self.reads.push(read.to_vec());
        self.n_reads += 1;
        Ok(())
    }

    pub fn consensus(&self) -> Result<Consensus, PoaError> {
        if self.n_reads < self.config.min_reads {
            return Err(PoaError::InsufficientDepth {
                got: self.n_reads,
                min: self.config.min_reads,
            });
        }

        // With only the seed, every edge has weight 1 so (weight-1)=0 everywhere.
        // heaviest_path would return the single first node. Return topology directly.
        if self.n_reads == 1 {
            let (topo, _) = topological_order(&self.nodes);
            let sequence: Vec<u8> = topo.iter().map(|&idx| self.nodes[idx].base).collect();
            let coverage: Vec<u32> = topo.iter().map(|_| 1).collect();
            let graph_stats = compute_stats(&self.nodes, self.config.min_allele_freq, self.n_reads);
            return Ok(Consensus {
                sequence,
                coverage,
                graph_stats,
            });
        }

        let (topo, rank_of) = topological_order(&self.nodes);
        let min_cov = self.min_coverage();

        let filtered: Vec<(usize, u8)> = match self.config.consensus_mode {
            ConsensusMode::HeaviestPath => {
                let path = heaviest_path(&self.nodes, &topo, &rank_of);

                // Trim leading nodes below min_cov
                let start = path
                    .iter()
                    .position(|&(node_idx, _)| self.nodes[node_idx].coverage >= min_cov)
                    .unwrap_or(0);

                // Trim trailing nodes below min_cov
                let end = path
                    .iter()
                    .rposition(|&(node_idx, _)| self.nodes[node_idx].coverage >= min_cov)
                    .map(|i| i + 1)
                    .unwrap_or(path.len());

                let effective = if start < end {
                    &path[start..end]
                } else {
                    &path[..]
                };

                // Filter interior minority detours
                effective
                    .iter()
                    .filter(|&&(node_idx, _)| self.nodes[node_idx].coverage >= min_cov)
                    .copied()
                    .collect()
            }
            ConsensusMode::MajorityFrequency => majority_frequency(&self.nodes, &topo, min_cov),
        };

        let sequence: Vec<u8> = filtered.iter().map(|&(_, base)| base).collect();
        let coverage: Vec<u32> = filtered
            .iter()
            .map(|&(node_idx, _)| self.nodes[node_idx].coverage)
            .collect();

        let graph_stats = compute_stats(&self.nodes, self.config.min_allele_freq, self.n_reads);

        Ok(Consensus {
            sequence,
            coverage,
            graph_stats,
        })
    }

    pub fn stats(&self) -> GraphStats {
        compute_stats(&self.nodes, self.config.min_allele_freq, self.n_reads)
    }

    /// Number of long-unbanded warnings emitted during `add_read` calls.
    /// Non-zero means at least one read exceeded 1 kb with `band_width = 0`.
    pub fn warnings_emitted(&self) -> usize {
        self.warnings
    }

    /// Build one consensus per detected allele by partitioning reads at the
    /// strongest bubble and running a fresh POA for each allele group.
    ///
    /// Returns a single-element `Vec` when no bubble is found (homozygous site
    /// or insufficient allele frequency). Returns `Err(InsufficientDepth)` when
    /// any allele group has fewer reads than `config.min_reads`.
    pub fn consensus_multi(&self) -> Result<Vec<Consensus>, PoaError> {
        if self.n_reads < self.config.min_reads {
            return Err(PoaError::InsufficientDepth {
                got: self.n_reads,
                min: self.config.min_reads,
            });
        }

        let (topo, _) = topological_order(&self.nodes);
        let bubbles = find_bubbles(
            &self.nodes,
            &topo,
            self.n_reads,
            self.config.min_allele_freq,
        );

        if bubbles.is_empty() {
            return Ok(vec![self.consensus()?]);
        }

        // Pick the bubble whose minor arm has the most support (strongest het signal).
        let (entry, arm_starts) = bubbles
            .iter()
            .max_by_key(|(entry, arms)| {
                arms.iter()
                    .filter_map(|&arm| self.edge_reads.get(&(*entry, arm)))
                    .map(|v| v.len())
                    .min()
                    .unwrap_or(0)
            })
            .unwrap();

        let groups = partition_reads_by_bubble(&self.edge_reads, *entry, arm_starts, self.n_reads);

        if groups.len() < 2 {
            return Ok(vec![self.consensus()?]);
        }

        let mut results = Vec::with_capacity(groups.len());
        for group in &groups {
            if group.len() < self.config.min_reads {
                return Err(PoaError::InsufficientDepth {
                    got: group.len(),
                    min: self.config.min_reads,
                });
            }
            let seed_slot = choose_seed(group, &self.reads);
            let seed = &self.reads[group[seed_slot]];
            let mut sub = PoaGraph::new(seed, self.config.clone())?;
            for (slot, &read_idx) in group.iter().enumerate() {
                if slot == seed_slot {
                    continue;
                }
                sub.add_read(&self.reads[read_idx])?;
            }
            results.push(sub.consensus()?);
        }

        Ok(results)
    }

    fn min_coverage(&self) -> u32 {
        if self.config.min_coverage_fraction > 0.0 {
            ((self.n_reads as f64 * self.config.min_coverage_fraction).ceil() as u32).max(1)
        } else if self.n_reads <= 1 {
            1
        } else {
            // Majority (n/2 + 1), minimum 2
            ((self.n_reads / 2 + 1).max(2)) as u32
        }
    }
}
