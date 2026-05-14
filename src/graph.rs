use crate::config::{AlignmentMode, ConsensusMode, PoaConfig};
use crate::error::PoaError;
use crate::types::{Consensus, CoverageGap, GapKind, GraphStats};
use std::collections::HashMap;

#[cfg(test)]
use std::cell::Cell as StdCell;

// Per-thread diagonal-skip counters, only compiled in test builds.
// SKIP_COUNTER: total diagonal skips fired; NODE_COUNTER: total non-source nodes.
#[cfg(test)]
thread_local! {
    pub(crate) static SKIP_COUNTER: StdCell<u64> = StdCell::new(0);
    pub(crate) static NODE_COUNTER: StdCell<u64> = StdCell::new(0);
}

#[cfg(test)]
pub(crate) fn reset_skip_counters() {
    SKIP_COUNTER.with(|c| c.set(0));
    NODE_COUNTER.with(|c| c.set(0));
}

#[cfg(test)]
pub(crate) fn skip_rate() -> f64 {
    let skips = SKIP_COUNTER.with(|c| c.get());
    let nodes = NODE_COUNTER.with(|c| c.get());
    if nodes == 0 {
        0.0
    } else {
        skips as f64 / nodes as f64
    }
}

// ─── Sentinel ────────────────────────────────────────────────────────────────

/// "Not yet filled" sentinel.
const UNSET: i32 = i32::MIN / 2;

#[inline]
fn safe_add(a: i32, b: i32) -> i32 {
    if a == UNSET {
        UNSET
    } else {
        a.saturating_add(b)
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
    /// Cached heaviest-path spine, recomputed adaptively rather than every read.
    cached_spine: Vec<(usize, u8, i32)>,
    /// `n_reads` at the time the spine was last recomputed.
    spine_updated_at: usize,
    /// Current recompute interval; doubles when the spine is stable, resets
    /// to 1 when it changes significantly.
    spine_interval: usize,
}

// ─── DP cell ─────────────────────────────────────────────────────────────────

/// Sentinel topo index meaning "came from the virtual start node".
const VIRTUAL: usize = usize::MAX;

#[derive(Clone, Copy, PartialEq, Eq)]
enum State {
    M,
    I,
    D,
}

#[derive(Clone, Copy)]
struct Cell {
    score: i32,
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

/// Measure the length of an arm starting at `start` by walking single-successor chains.
///
/// Stops at `max_depth`, a fork (multiple out-edges), or a reconvergence point
/// (a node reached from multiple in-edges, after the first step).
/// The reconvergence node itself is NOT counted.
fn materialize_arm_len(nodes: &[Node], start: usize, max_depth: usize) -> usize {
    let mut len = 0;
    let mut cur = start;
    for _ in 0..max_depth {
        // After the first step, a node with multiple in-edges is a reconvergence.
        if len >= 1 && nodes[cur].in_edges.len() > 1 {
            break;
        }
        len += 1;
        match nodes[cur].out_edges.as_slice() {
            [(next, _)] => cur = *next,
            _ => break,
        }
    }
    len
}

// ─── Bubble detection ────────────────────────────────────────────────────────

/// Find the topo rank of the exit node for a bubble starting at `entry_t`.
///
/// Uses the path-budget algorithm: budget starts at `out_degree(entry)` — one open
/// path per arm.  For each subsequent node (in topo order) that is reachable from
/// entry within the bubble, budget decrements by the number of in-edges arriving
/// from already-seen bubble nodes, then increments by the node's own out-degree.
/// When budget reaches 0 all paths have reconverged — that node is the exit.
///
/// Returns `None` when the graph ends before reconvergence (open-ended bubble).
fn find_bubble_exit(nodes: &[Node], topo: &[usize], entry_t: usize) -> Option<usize> {
    let entry_node = topo[entry_t];
    let nn = nodes.len();
    let mut in_bubble = vec![false; nn];
    in_bubble[entry_node] = true;

    let mut outstanding = nodes[entry_node].out_edges.len();

    for t in (entry_t + 1)..topo.len() {
        let node_idx = topo[t];
        let bubble_in = nodes[node_idx]
            .in_edges
            .iter()
            .filter(|&&p| in_bubble[p])
            .count();

        if bubble_in == 0 {
            continue;
        }

        in_bubble[node_idx] = true;
        outstanding -= bubble_in;

        if outstanding == 0 {
            return Some(t);
        }

        outstanding += nodes[node_idx].out_edges.len();
    }

    None
}

/// For each topo rank, compute the bubble it belongs to as `Some((entry_t, exit_t))`,
/// or `None` for spine nodes outside any bubble.
///
/// Bubble nodes restrict the DP to a window sized by the bubble span; spine nodes
/// use a narrow diagonal window.  All nodes in `entry_t..=exit_t` share the same
/// entry's j position as their window anchor.
fn compute_bubble_ranges(nodes: &[Node], topo: &[usize]) -> Vec<Option<(usize, usize)>> {
    let n = topo.len();
    let mut ranges: Vec<Option<(usize, usize)>> = vec![None; n];
    let mut t = 0;
    while t < n {
        if nodes[topo[t]].out_edges.len() >= 2 {
            if let Some(exit_t) = find_bubble_exit(nodes, topo, t) {
                for bt in t..=exit_t {
                    ranges[bt] = Some((t, exit_t));
                }
                t = exit_t + 1;
            } else {
                // Open-ended bubble: mark everything remaining conservatively.
                for bt in t..n {
                    ranges[bt] = Some((t, n - 1));
                }
                break;
            }
        } else {
            t += 1;
        }
    }
    ranges
}

// ─── Stale-spine helpers ─────────────────────────────────────────────────────

/// Base-level differences between two spines: length delta + per-position
/// mismatches over the shared prefix.  Returns `usize::MAX` when `old` is
/// empty (forces the first recompute to be treated as fully unstable).
fn spine_diff(old: &[(usize, u8, i32)], new: &[(usize, u8, i32)]) -> usize {
    if old.is_empty() {
        return usize::MAX;
    }
    let len_diff = old.len().abs_diff(new.len());
    let base_diffs = old[..old.len().min(new.len())]
        .iter()
        .zip(new.iter())
        .filter(|(o, n)| o.1 != n.1)
        .count();
    base_diffs + len_diff
}

/// Spine is "stable" if it changed by this many bases or fewer.
const SPINE_STABLE_THRESHOLD: usize = 3;
/// Maximum recompute interval (reads between spine refreshes).
const SPINE_MAX_INTERVAL: usize = 32;

// ─── Lookahead resolve ───────────────────────────────────────────────────────

/// Query bases scored per arm when resolving a bubble entry.
/// Set to one repeat unit (5 bases covers AAGGG, GGGGCC, CAG etc.).
/// Larger windows increase false-positive risk at typical ONT error rates
/// because a single read error in the window can flip the decision.
const LOOKAHEAD_K: usize = 5;
/// Winning arm must beat the runner-up by at least this many score points.
/// One match+mismatch swing = 2 points, so MARGIN=2 requires one net
/// advantage position after accounting for any tie positions.
const LOOKAHEAD_MARGIN: i32 = 2;
/// Maximum nodes walked when materialising a single bubble arm.
const LOOKAHEAD_MAX_ARM_DEPTH: usize = 512;

/// Walks one bubble arm from `start_idx` forward, collecting node indices
/// until the bubble exit (topo rank `exit_t`) is reached.  The exit node
/// itself is NOT included — it runs its own windowed DP normally.
///
/// Returns `None` if the arm has internal branching (a nested bubble) or
/// exceeds the depth cap; both are signals to fall back to windowed DP.
fn collect_arm_nodes(
    nodes: &[Node],
    rank_of: &[usize],
    start_idx: usize,
    exit_t: usize,
) -> Option<Vec<usize>> {
    let mut path = Vec::new();
    let mut cur = start_idx;
    for _ in 0..LOOKAHEAD_MAX_ARM_DEPTH {
        if rank_of[cur] == exit_t {
            return Some(path);
        }
        path.push(cur);
        match nodes[cur].out_edges.as_slice() {
            [(next, _)] => cur = *next,
            _ => return None, // internal branch or dead end before exit
        }
    }
    None // depth cap exceeded
}

/// Scores the first `min(LOOKAHEAD_K, arm.len(), query[j..].len())` bases of
/// `arm` against `query` starting at position `j`.
fn score_arm_prefix(
    arm: &[usize],
    nodes: &[Node],
    query: &[u8],
    j: usize,
    match_score: i32,
    mismatch_score: i32,
) -> i32 {
    let len = arm
        .len()
        .min(LOOKAHEAD_K)
        .min(query.len().saturating_sub(j));
    arm[..len]
        .iter()
        .enumerate()
        .map(|(i, &nidx)| {
            if nodes[nidx].base == query[j + i] {
                match_score
            } else {
                mismatch_score
            }
        })
        .sum()
}

// ─── Bubble-aware DP alignment ────────────────────────────────────────────────

/// Half-width of the DP band applied to spine (non-bubble) nodes.
/// Covers ±SPINE_MARGIN query positions around the expected diagonal.
/// 50 accommodates up to ~10 % indel rate over a 500-node spine segment.
const SPINE_MARGIN: usize = 50;

/// Align `query` against the graph using spine-guided affine-gap DAG DP.
///
/// **Diagonal skip** — O(1) fast path for consecutive exact matches along the
/// heaviest-path spine.  Fires whenever the spine node's base matches the query
/// base AND the M score at the predecessor dominates any open I/D state.  After
/// each read is added the spine converges toward the true consensus, so later
/// reads align almost entirely via diagonal skip with only small local-DP
/// bursts at divergence points.
///
/// **j-range restriction** — spine nodes get a narrow ±[`SPINE_MARGIN`] window;
/// bubble nodes get a window anchored on the bubble entry's j position extended
/// by the bubble span.  Cells outside the window remain UNSET and are invisible
/// to traceback.
fn align(
    nodes: &[Node],
    topo: &[usize],
    rank_of: &[usize],
    spine: &[(usize, u8, i32)],
    query: &[u8],
    cfg: &PoaConfig,
) -> Result<Vec<AlignOp>, PoaError> {
    let n = topo.len();
    let l = query.len();
    let go = cfg.gap_open;
    let ge = cfg.gap_extend;
    let semi = cfg.alignment_mode == AlignmentMode::SemiGlobal;

    // Flat DP tables: index = t * (l+1) + j
    let size = n * (l + 1);
    let mut m = vec![Cell::unset(); size];
    let mut ins = vec![Cell::unset(); size];
    let mut del = vec![Cell::unset(); size];

    // Spine membership: for each node index, the previous spine node (if on spine).
    // Used to gate the relaxed diagonal skip.
    let nn = nodes.len();
    let mut on_spine = vec![false; nn];
    let mut spine_prev: Vec<Option<usize>> = vec![None; nn];
    for (sp, &(node_idx, _, _)) in spine.iter().enumerate() {
        on_spine[node_idx] = true;
        if sp > 0 {
            spine_prev[node_idx] = Some(spine[sp - 1].0);
        }
    }

    // Bubble membership: Some((entry_t, exit_t)) or None for spine.
    let bubble_ranges = compute_bubble_ranges(nodes, topo);
    // j position of each bubble's entry node's predecessors — filled dynamically.
    let mut bubble_entry_j = vec![0usize; n];

    // Track best_j per row for diagonal skip.
    let mut best_j_per_t = vec![0usize; n];

    // Nodes marked here are on a losing bubble arm identified by lookahead
    // resolve.  Their entire windowed DP is skipped; their cells stay UNSET so
    // the exit node naturally reads only from the winning arm's predecessors.
    let mut lookahead_skip = vec![false; nn];

    #[inline]
    fn idx(t: usize, j: usize, l: usize) -> usize {
        t * (l + 1) + j
    }

    for (t, &node_idx) in topo.iter().enumerate() {
        // Lookahead: skip nodes on losing bubble arms entirely.
        if lookahead_skip[node_idx] {
            continue;
        }

        let node_base = nodes[node_idx].base;
        let is_source = nodes[node_idx].in_edges.is_empty();

        // ── Diagonal skip ──────────────────────────────────────────────────────
        // O(1) fast path for consecutive exact matches along the spine.
        //
        // Correct formulation: predecessor has consumed `bj` query bases
        // (M[pt][bj] is its best cell).  We advance by one — consuming query[bj]
        // (0-indexed) — and record M[t][bj+1] = M[pt][bj] + match_score.
        // best_j[t] = bj+1, so the NEXT spine node will find M[t][bj+1] set and
        // fire its own skip.  This allows runs of consecutive diagonal skips.
        //
        // Guards:
        //   - node must be on the spine (so we know which predecessor to trust)
        //   - exactly one in-edge, from the previous spine node (no competing arm)
        //   - exactly one out-edge (not a bubble entry — those need full DP so the
        //     non-spine arms can read correct predecessor scores)
        //   - no open insert/delete state at predecessor's best column
        #[cfg(test)]
        if !is_source {
            NODE_COUNTER.with(|c| c.set(c.get() + 1));
        }

        if !is_source && on_spine[node_idx] {
            if let Some(prev_sp) = spine_prev[node_idx] {
                if nodes[node_idx].in_edges == [prev_sp] && nodes[node_idx].out_edges.len() == 1 {
                    let pt = rank_of[prev_sp];
                    let bj = best_j_per_t[pt];
                    if bj < l && node_base == query[bj] {
                        let m_prev = m[idx(pt, bj, l)].score;
                        let i_prev = ins[idx(pt, bj, l)].score;
                        let d_prev = del[idx(pt, bj, l)].score;
                        // The source node always has I[source][j] set from global
                        // alignment initialization (leading-insert penalties).
                        // It's safe to skip through source if M clearly dominates —
                        // I[source][bj] is just the penalty for a leading insertion,
                        // not an ongoing gap run.  For all non-source predecessors
                        // we require I/D to be strictly UNSET: an ongoing gap run
                        // must not be cut off by the skip.
                        let pred_is_source = nodes[topo[pt]].in_edges.is_empty();
                        let i_ok = i_prev == UNSET || (pred_is_source && m_prev > i_prev);
                        let d_ok = d_prev == UNSET || (pred_is_source && m_prev > d_prev);
                        if m_prev != UNSET && i_ok && d_ok {
                            let score = m_prev + cfg.match_score;
                            m[idx(t, bj + 1, l)] = Cell { score, pred_t: pt };
                            best_j_per_t[t] = bj + 1;
                            #[cfg(test)]
                            SKIP_COUNTER.with(|c| c.set(c.get() + 1));
                            continue;
                        }
                    }
                }
            }
        }

        // ── j = 0: delete-only ──────────────────────────────────────────────────
        {
            let ix0 = idx(t, 0, l);
            let (mut best, mut best_pred) = (UNSET, 0usize);
            if is_source {
                let val = go + ge;
                if val > best {
                    best = val;
                    best_pred = VIRTUAL;
                }
            }
            for &p in &nodes[node_idx].in_edges {
                let pt = rank_of[p];
                let vm = safe_add(m[idx(pt, 0, l)].score, go + ge);
                if vm > best {
                    best = vm;
                    best_pred = pt;
                }
                let vd = safe_add(del[idx(pt, 0, l)].score, ge);
                if vd > best {
                    best = vd;
                    best_pred = pt;
                }
            }
            if best != UNSET {
                del[ix0] = Cell {
                    score: best,
                    pred_t: best_pred,
                };
            }
        }

        // ── j range: narrow window for both spine and bubble nodes ────────────
        // The window is centred on the expected diagonal position.
        //   Spine node  : j_center = predecessor's best_j + 1  (±SPINE_MARGIN)
        //   Bubble node : j_center = bubble-entry's predecessor j
        //                 hi  extended by bubble span so all arm lengths fit
        // Source nodes always use the full range.
        // Cells outside the range remain UNSET and traceback skips them.
        let j_pred_max = if is_source {
            0
        } else {
            nodes[node_idx]
                .in_edges
                .iter()
                .map(|&p| best_j_per_t[rank_of[p]])
                .max()
                .unwrap_or(0)
        };

        let (j_lo, j_hi) = if is_source {
            (1usize, l)
        } else {
            match bubble_ranges[t] {
                Some((entry_t, exit_t)) => {
                    // Record this bubble's entry-predecessor j the first time we see it.
                    if t == entry_t {
                        bubble_entry_j[entry_t] = j_pred_max;
                    }
                    let bej = bubble_entry_j[entry_t];
                    let bubble_span = exit_t.saturating_sub(entry_t) + 1;
                    let lo = bej.saturating_sub(SPINE_MARGIN).max(1);
                    let hi = (bej + bubble_span + SPINE_MARGIN).min(l);
                    (lo, hi)
                }
                None => {
                    // Spine: narrow diagonal window.
                    let j_center = j_pred_max.saturating_add(1);
                    let lo = j_center.saturating_sub(SPINE_MARGIN).max(1);
                    let hi = (j_center + SPINE_MARGIN).min(l);
                    (lo, hi)
                }
            }
        };

        // ── j = j_lo..=j_hi ───────────────────────────────────────────────────
        let mut row_best_j = 0usize;
        let mut row_best_score = UNSET;

        for j in j_lo..=j_hi {
            let q_base = query[j - 1];
            let sc = if node_base == q_base {
                cfg.match_score
            } else {
                cfg.mismatch_score
            };
            let ixcur = idx(t, j, l);

            // M[t][j]
            {
                let (mut best, mut best_pred) = (UNSET, VIRTUAL);
                // Free-start at j=1: applies when j_lo==1 (source or bubble nodes).
                if j == 1 && (is_source || semi) {
                    if sc > best {
                        best = sc;
                        best_pred = VIRTUAL;
                    }
                }
                if is_source && j > 1 {
                    let val = safe_add(go + (j as i32 - 1) * ge, sc);
                    if val != UNSET && val > best {
                        best = val;
                        best_pred = VIRTUAL;
                    }
                }
                for &p in &nodes[node_idx].in_edges {
                    let pt = rank_of[p];
                    let vm = safe_add(m[idx(pt, j - 1, l)].score, sc);
                    if vm != UNSET && vm > best {
                        best = vm;
                        best_pred = pt;
                    }
                    let vi = safe_add(ins[idx(pt, j - 1, l)].score, sc);
                    if vi != UNSET && vi > best {
                        best = vi;
                        best_pred = pt;
                    }
                    let vd = safe_add(del[idx(pt, j - 1, l)].score, sc);
                    if vd != UNSET && vd > best {
                        best = vd;
                        best_pred = pt;
                    }
                }
                if best != UNSET {
                    m[ixcur] = Cell {
                        score: best,
                        pred_t: best_pred,
                    };
                    if best > row_best_score {
                        row_best_score = best;
                        row_best_j = j;
                    }
                }
            }

            // I[t][j]
            {
                let ixprev = idx(t, j - 1, l);
                let (mut best, mut best_pred) = (UNSET, VIRTUAL);
                if is_source && j == 1 {
                    let val = go + ge;
                    if val > best {
                        best = val;
                        best_pred = VIRTUAL;
                    }
                }
                let vm = safe_add(m[ixprev].score, go + ge);
                if vm != UNSET && vm > best {
                    best = vm;
                    best_pred = t;
                }
                let vi = safe_add(ins[ixprev].score, ge);
                if vi != UNSET && vi > best {
                    best = vi;
                    best_pred = t;
                }
                if best != UNSET {
                    ins[ixcur] = Cell {
                        score: best,
                        pred_t: best_pred,
                    };
                }
            }

            // D[t][j]
            {
                let (mut best, mut best_pred) = (UNSET, 0usize);
                for &p in &nodes[node_idx].in_edges {
                    let pt = rank_of[p];
                    let vm = safe_add(m[idx(pt, j, l)].score, go + ge);
                    if vm != UNSET && vm > best {
                        best = vm;
                        best_pred = pt;
                    }
                    let vi = safe_add(ins[idx(pt, j, l)].score, go + ge);
                    if vi != UNSET && vi > best {
                        best = vi;
                        best_pred = pt;
                    }
                    let vd = safe_add(del[idx(pt, j, l)].score, ge);
                    if vd != UNSET && vd > best {
                        best = vd;
                        best_pred = pt;
                    }
                }
                if best != UNSET {
                    del[ixcur] = Cell {
                        score: best,
                        pred_t: best_pred,
                    };
                }
            }
        }

        if row_best_score != UNSET {
            best_j_per_t[t] = row_best_j;
        } else {
            // No M cell set: estimate j from predecessor + 1 (diagonal advance).
            best_j_per_t[t] = if is_source {
                1
            } else {
                j_pred_max.saturating_add(1).min(l)
            };
        }

        // ── Lookahead resolve ─────────────────────────────────────────────────
        // At a bubble entry node: score each arm's first LOOKAHEAD_K bases
        // against the query.  If one arm beats the runner-up by LOOKAHEAD_MARGIN
        // points, mark every losing arm's nodes in `lookahead_skip` so their
        // windowed DP is skipped.  The winning arm runs its normal windowed DP
        // (the existing j-window is already correctly anchored to bubble_entry_j).
        // The exit node then reads from the winning arm's cells — losing arms
        // are all UNSET and invisible to traceback.
        if let Some((entry_t, exit_t)) = bubble_ranges[t] {
            if t == entry_t && nodes[node_idx].out_edges.len() >= 2 {
                let j_entry = best_j_per_t[t];
                if j_entry < l {
                    // Materialise each arm; bail if any is internally branched.
                    let mut all_arms: Vec<Vec<usize>> = Vec::new();
                    let mut complex = false;
                    for &(arm_start, _) in &nodes[node_idx].out_edges {
                        match collect_arm_nodes(nodes, rank_of, arm_start, exit_t) {
                            Some(arm) => all_arms.push(arm),
                            None => {
                                complex = true;
                                break;
                            }
                        }
                    }

                    // Only fire lookahead when the longest arm provides the
                    // full LOOKAHEAD_K bases to score.  Short arms (< K nodes)
                    // are unreliable: a single read error in the window can flip
                    // the decision, causing incorrect arm elimination.  Requiring
                    // the full K bases means only genuine structural bubbles
                    // (≥ one repeat unit long) trigger lookahead; error-induced
                    // single-node arms always fall back to windowed DP.
                    let max_scored = all_arms
                        .iter()
                        .map(|arm| {
                            arm.len()
                                .min(LOOKAHEAD_K)
                                .min(query.len().saturating_sub(j_entry))
                        })
                        .max()
                        .unwrap_or(0);

                    if !complex && all_arms.len() >= 2 && max_scored >= LOOKAHEAD_K {
                        let scores: Vec<i32> = all_arms
                            .iter()
                            .map(|arm| {
                                score_arm_prefix(
                                    arm,
                                    nodes,
                                    query,
                                    j_entry,
                                    cfg.match_score,
                                    cfg.mismatch_score,
                                )
                            })
                            .collect();

                        let best = *scores.iter().max().unwrap();
                        let winner_count = scores.iter().filter(|&&s| s == best).count();

                        if winner_count == 1 {
                            let second_best = scores
                                .iter()
                                .copied()
                                .filter(|&s| s != best)
                                .max()
                                .unwrap_or(i32::MIN);

                            if best - second_best >= LOOKAHEAD_MARGIN {
                                let winner = scores.iter().position(|&s| s == best).unwrap();
                                for (arm_idx, arm) in all_arms.iter().enumerate() {
                                    if arm_idx != winner {
                                        for &losing_node in arm {
                                            lookahead_skip[losing_node] = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // ── Find best terminal cell at column l ─────────────────────────────────────
    let terminal_best = (0..n)
        .flat_map(|t| {
            let sm = m[idx(t, l, l)].score;
            let si = ins[idx(t, l, l)].score;
            let sd = del[idx(t, l, l)].score;
            let best_sc = [sm, si, sd].into_iter().filter(|&s| s != UNSET).max();
            best_sc.map(|sc| {
                let st = if sm == sc {
                    State::M
                } else if si == sc {
                    State::I
                } else {
                    State::D
                };
                (t, st, sc)
            })
        })
        .max_by_key(|&(_, _, sc)| sc)
        .map(|(t, s, _)| (t, s));

    let (best_t, best_state) = match terminal_best {
        Some(result) => result,
        None => (0, State::M),
    };

    // ── Traceback ─────────────────────────────────────────────────────────────
    let mut ops: Vec<AlignOp> = Vec::with_capacity(l + n / 4);
    let mut t = best_t;
    let mut j = l;
    let mut cur_state = best_state;

    loop {
        let cell = match cur_state {
            State::M => m[idx(t, j, l)],
            State::I => ins[idx(t, j, l)],
            State::D => del[idx(t, j, l)],
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
                cur_state = best_prev_state(&m, &ins, &del, t, j, l);
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
                cur_state = if m[idx(t, j, l)].score >= ins[idx(t, j, l)].score {
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
                cur_state = best_prev_state(&m, &ins, &del, t, j, l);
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
fn best_prev_state(m: &[Cell], ins: &[Cell], del: &[Cell], t: usize, j: usize, l: usize) -> State {
    let sm = m[t * (l + 1) + j].score;
    let si = ins[t * (l + 1) + j].score;
    let sd = del[t * (l + 1) + j].score;
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

fn heaviest_path(nodes: &[Node], topo: &[usize], rank_of: &[usize]) -> Vec<(usize, u8, i32)> {
    let n = topo.len();
    let mut cum: Vec<(i64, Option<usize>, i32)> = vec![(0, None, 0); n];

    for t in 0..n {
        let node_idx = topo[t];
        let curr = cum[t].0;
        for &(succ_idx, weight) in &nodes[node_idx].out_edges {
            let succ_t = rank_of[succ_idx];
            let candidate = curr + (weight - 1) as i64;
            if candidate > cum[succ_t].0 {
                cum[succ_t] = (candidate, Some(t), weight);
            }
        }
    }

    let max_cum = (0..n).map(|t| cum[t].0).max().unwrap_or(0);
    let best_t = (0..n).find(|&t| cum[t].0 == max_cum).unwrap_or(0);

    let mut path: Vec<(usize, u8, i32)> = Vec::new();
    let mut t = best_t;
    loop {
        let node_idx = topo[t];
        let w = if cum[t].1.is_none() {
            nodes[node_idx].coverage as i32
        } else {
            cum[t].2
        };
        path.push((node_idx, nodes[node_idx].base, w));
        match cum[t].1 {
            None => break,
            Some(pred_t) => t = pred_t,
        }
    }
    path.reverse();
    path
}

// ─── Majority-frequency consensus ────────────────────────────────────────────

fn majority_frequency(nodes: &[Node], topo: &[usize], min_cov: u32) -> Vec<(usize, u8, i32)> {
    topo.iter()
        .copied()
        .filter(|&idx| {
            let cov = nodes[idx].coverage;
            let del = nodes[idx].delete_count;
            let total = cov + del;
            total >= min_cov && cov * 2 >= total
        })
        .map(|idx| (idx, nodes[idx].base, nodes[idx].coverage as i32))
        .collect()
}

// ─── Graph statistics ─────────────────────────────────────────────────────────

fn compute_stats(nodes: &[Node], min_allele_freq: f64, n_reads: usize) -> GraphStats {
    let node_count = nodes.len();
    let edge_count: usize = nodes.iter().map(|nd| nd.out_edges.len()).sum();

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

    let single_support = nodes.iter().filter(|nd| nd.coverage == 1).count();
    let single_support_fraction = if node_count == 0 {
        0.0
    } else {
        single_support as f64 / node_count as f64
    };

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
            qualifying.sort_unstable_by(|a, b| b.cmp(a));
            max_bubble_depth = max_bubble_depth.max(qualifying[1] as usize);
        }
    }

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

// ─── Coverage gap detection ───────────────────────────────────────────────────

fn detect_coverage_gaps(coverage: &[u32]) -> Vec<CoverageGap> {
    let first = coverage.iter().position(|&c| c >= 2);
    let last = coverage.iter().rposition(|&c| c >= 2);
    let (first, last) = match (first, last) {
        (Some(f), Some(l)) if f < l => (f, l),
        _ => return vec![],
    };
    let mut gaps = Vec::new();
    let mut gap_start: Option<usize> = None;
    for (offset, &cov) in coverage[(first + 1)..last].iter().enumerate() {
        let i = first + 1 + offset;
        if cov < 2 {
            gap_start.get_or_insert(i);
        } else if let Some(s) = gap_start.take() {
            gaps.push(CoverageGap {
                start: s,
                end: i,
                kind: GapKind::Spanning,
            });
        }
    }
    if let Some(s) = gap_start {
        gaps.push(CoverageGap {
            start: s,
            end: last,
            kind: GapKind::Spanning,
        });
    }
    gaps
}

// ─── Multi-allele helpers ─────────────────────────────────────────────────────

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

fn partition_reads_by_bubble(
    edge_reads: &HashMap<(usize, usize), Vec<u32>>,
    entry_node: usize,
    arm_starts: &[usize],
    n_reads: usize,
) -> Vec<Vec<usize>> {
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
            cached_spine: Vec::new(),
            spine_updated_at: 0,
            spine_interval: 1,
        })
    }

    pub fn add_read(&mut self, read: &[u8]) -> Result<(), PoaError> {
        if read.is_empty() {
            return Err(PoaError::EmptyInput);
        }

        let (topo, rank_of) = topological_order(&self.nodes);

        // Refresh the spine when the cache is empty or the interval has elapsed.
        // The interval doubles each time the spine is stable (≤ SPINE_STABLE_THRESHOLD
        // base changes), and resets to 1 when it changes significantly.  The final
        // consensus always recomputes the heaviest path from scratch, so a stale
        // spine only affects alignment speed, never correctness.
        let reads_since_update = self.n_reads.saturating_sub(self.spine_updated_at);
        if self.cached_spine.is_empty() || reads_since_update >= self.spine_interval {
            let new_spine = heaviest_path(&self.nodes, &topo, &rank_of);
            let diff = spine_diff(&self.cached_spine, &new_spine);
            self.cached_spine = new_spine;
            self.spine_updated_at = self.n_reads;
            if diff <= SPINE_STABLE_THRESHOLD {
                self.spine_interval = (self.spine_interval * 2).min(SPINE_MAX_INTERVAL);
            } else {
                self.spine_interval = 1;
            }
        }

        let ops = align(
            &self.nodes,
            &topo,
            &rank_of,
            &self.cached_spine,
            read,
            &self.config,
        )?;
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

        if self.n_reads == 1 {
            let (topo, _) = topological_order(&self.nodes);
            let sequence: Vec<u8> = topo.iter().map(|&idx| self.nodes[idx].base).collect();
            let coverage: Vec<u32> = topo.iter().map(|_| 1).collect();
            let path_weights: Vec<i32> = topo.iter().map(|_| 1).collect();
            let graph_stats = compute_stats(&self.nodes, self.config.min_allele_freq, self.n_reads);
            return Ok(Consensus {
                sequence,
                coverage,
                path_weights,
                n_reads: 1,
                graph_stats,
                gaps: vec![],
            });
        }

        let (topo, rank_of) = topological_order(&self.nodes);
        let min_cov = self.min_coverage();

        let filtered: Vec<(usize, u8, i32)> = match self.config.consensus_mode {
            ConsensusMode::HeaviestPath => {
                let path = heaviest_path(&self.nodes, &topo, &rank_of);

                let start = path
                    .iter()
                    .position(|&(node_idx, _, _)| self.nodes[node_idx].coverage >= min_cov)
                    .unwrap_or(0);

                let end = path
                    .iter()
                    .rposition(|&(node_idx, _, _)| self.nodes[node_idx].coverage >= min_cov)
                    .map(|i| i + 1)
                    .unwrap_or(path.len());

                let effective = if start < end {
                    &path[start..end]
                } else {
                    &path[..]
                };

                effective
                    .iter()
                    .filter(|&&(node_idx, _, _)| self.nodes[node_idx].coverage >= min_cov)
                    .copied()
                    .collect()
            }
            ConsensusMode::MajorityFrequency => majority_frequency(&self.nodes, &topo, min_cov),
        };

        let sequence: Vec<u8> = filtered.iter().map(|&(_, base, _)| base).collect();
        let coverage: Vec<u32> = filtered
            .iter()
            .map(|&(node_idx, _, _)| self.nodes[node_idx].coverage)
            .collect();
        let path_weights: Vec<i32> = filtered.iter().map(|&(_, _, w)| w).collect();

        let graph_stats = compute_stats(&self.nodes, self.config.min_allele_freq, self.n_reads);
        let gaps = detect_coverage_gaps(&coverage);

        Ok(Consensus {
            sequence,
            coverage,
            path_weights,
            n_reads: self.n_reads,
            graph_stats,
            gaps,
        })
    }

    pub fn stats(&self) -> GraphStats {
        compute_stats(&self.nodes, self.config.min_allele_freq, self.n_reads)
    }

    /// Number of long-unbanded warnings emitted during `add_read` calls.
    /// With the new bubble-aware aligner this always returns 0 (warnings removed).
    pub fn warnings_emitted(&self) -> usize {
        self.warnings
    }

    /// Return all edge weights in the graph (one entry per directed edge).
    pub fn edge_weights(&self) -> Vec<i32> {
        self.nodes
            .iter()
            .flat_map(|n| n.out_edges.iter().map(|&(_, w)| w))
            .collect()
    }

    /// Return per-node coverage in topological order.
    pub fn node_coverages(&self) -> Vec<u32> {
        let (topo, _) = topological_order(&self.nodes);
        topo.iter().map(|&i| self.nodes[i].coverage).collect()
    }

    /// Align `read` into the current graph and return the alignment operations
    /// without modifying the graph.
    pub fn align_read_ops(
        &self,
        read: &[u8],
    ) -> Result<(Vec<AlignOp>, usize, Vec<usize>), PoaError> {
        let (topo, rank_of) = topological_order(&self.nodes);
        let spine = heaviest_path(&self.nodes, &topo, &rank_of);
        let ops = align(&self.nodes, &topo, &rank_of, &spine, read, &self.config)?;
        Ok((ops, 0, rank_of))
    }

    /// Like [`align_read_ops`] but kept for compatibility.
    pub fn align_read_ops_unbanded(
        &self,
        read: &[u8],
    ) -> Result<(Vec<AlignOp>, Vec<usize>), PoaError> {
        let (topo, rank_of) = topological_order(&self.nodes);
        let spine = heaviest_path(&self.nodes, &topo, &rank_of);
        let ops = align(&self.nodes, &topo, &rank_of, &spine, read, &self.config)?;
        Ok((ops, rank_of))
    }

    /// Number of nodes currently in the graph.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// For each node, return `(out_edges_count, max_out_edge_weight, min_out_edge_weight)`.
    pub fn node_out_edge_info(&self) -> Vec<(usize, i32, i32)> {
        self.nodes
            .iter()
            .map(|n| {
                let count = n.out_edges.len();
                let max_w = n.out_edges.iter().map(|&(_, w)| w).max().unwrap_or(0);
                let min_w = n.out_edges.iter().map(|&(_, w)| w).min().unwrap_or(0);
                (count, max_w, min_w)
            })
            .collect()
    }

    /// Return arm lengths for each bubble node with 2+ qualifying arms (weight >= threshold).
    ///
    /// Arm length is measured by walking the single-successor chain from each arm start,
    /// stopping at reconvergence points (nodes with multiple in-edges after the first step).
    pub fn bubble_arm_lengths(
        &self,
        weight_threshold: i32,
        arm_len_threshold: usize,
    ) -> Vec<(usize, Vec<usize>)> {
        let (topo, _) = topological_order(&self.nodes);
        let mut result = Vec::new();
        for (t, &node_idx) in topo.iter().enumerate() {
            let qualifying: Vec<usize> = self.nodes[node_idx]
                .out_edges
                .iter()
                .filter(|&&(_, w)| w >= weight_threshold)
                .map(|&(succ, _)| succ)
                .collect();
            if qualifying.len() >= 2 {
                // Measure each arm by walking single-successor chains.
                let arm_lens: Vec<usize> = qualifying
                    .iter()
                    .map(|&arm_start| materialize_arm_len(&self.nodes, arm_start, 500))
                    .collect();
                let min_len = arm_lens.iter().copied().min().unwrap_or(0);
                if min_len >= arm_len_threshold {
                    result.push((t, arm_lens));
                }
            }
        }
        result
    }

    /// Build one consensus per detected allele.
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
            ((self.n_reads / 2 + 1).max(2)) as u32
        }
    }
}
