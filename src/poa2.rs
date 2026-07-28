//! Clean partial-order alignment consensus engine (abPOA-style rebuild).
//!
//! Deliberately mirrors abPOA's proven-clean core — full partial-order DP (one
//! DP row per graph node, every node relaxes ALL its in-edges), heaviest
//! bundling for consensus — with none of the accumulated machinery from
//! `graph.rs` (spine-guided band, diagonal-skip, slide-lock, bypass edges,
//! interior filter). The value the crate adds lives ABOVE this engine (analysis,
//! CI, diagnostics); the engine itself just aims to match abPOA quality.
//!
//! Increment 1 (this file): data model, seed, full-PO affine **global** DP,
//! read integration with SPOA/abPOA-style **node fusion** (match → reuse node;
//! mismatch → reuse an aligned sibling carrying that base, else create a new
//! substitution node and cross-link it into the alignment column; insert → new
//! node; delete → skip), heaviest-bundling consensus. Node fusion is in the core
//! deliberately: full DP does NOT obviate it (DP decides the alignment; the
//! graph-update must still represent a mismatch as a substitution node, or the
//! read's base is dropped and the seed's minority base is wrongly kept — see
//! `tests_mismatch::majority_substitution_over_seed_base`). abPOA/SPOA both do
//! this via aligned-node columns.
//!
//! Not yet (later increments, for scale/real data): adaptive static-diagonal-
//! union band + semi-global free ends (full global DP is O(nodes×readlen) and
//! fixed-end — fine for unit tests, too heavy/rigid for multi-kb partial reads);
//! then real-data validation vs abPOA; then two-piece convex gaps and SIMD.

use crate::config::{AlignmentMode, PoaConfig};
use std::collections::HashMap;

const NEG_INF: i32 = i32::MIN / 4; // /4 so additions can't overflow

/// Consensus coverage floor (mirrors the crate's `coverage_threshold`): a
/// consensus node below this at a boundary is trimmed. Default fraction 0.0
/// falls back to a strict majority `(n/2+1).max(2)`.
fn coverage_threshold(population: usize, min_coverage_fraction: f64) -> u32 {
    if min_coverage_fraction > 0.0 {
        ((population as f64 * min_coverage_fraction).ceil() as u32).max(1)
    } else if population <= 1 {
        1
    } else {
        ((population / 2 + 1).max(2)) as u32
    }
}

/// Read a banded matrix cell (M or I). Returns NEG_INF for j=0 or out-of-band.
#[inline]
fn bget(flat: &[i32], off: &[usize], lo: &[usize], hi: &[usize], r: usize, j: usize) -> i32 {
    if j == 0 || j < lo[r] || j > hi[r] {
        NEG_INF
    } else {
        flat[off[r] + (j - lo[r])]
    }
}

/// Read a banded D cell. j=0 lives in the separate `d0` column.
#[inline]
fn bgetd(
    flat: &[i32],
    d0: &[i32],
    off: &[usize],
    lo: &[usize],
    hi: &[usize],
    r: usize,
    j: usize,
) -> i32 {
    if j == 0 {
        d0[r]
    } else if j < lo[r] || j > hi[r] {
        NEG_INF
    } else {
        flat[off[r] + (j - lo[r])]
    }
}

/// An out-edge and how many reads traversed it (Match/Insert founding traffic).
#[derive(Clone, Copy)]
struct Edge {
    to: usize,
    weight: i32,
}

struct Node {
    base: u8,
    out: Vec<Edge>,
    inc: Vec<usize>,     // predecessor node indices
    aligned: Vec<usize>, // sibling nodes in the same MSA column (alternative bases)
    cov: u32,            // reads that Matched/founded this node (per-base depth)
    del: u32,            // reads that Deleted (skipped) this node
}

/// Alignment op from traceback (forward read order).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Op {
    Match(usize), // read base aligned to this node (match or mismatch)
    Ins(u8),      // read base with no graph node
    Del(usize),   // graph node skipped by the read
}

pub struct Poa {
    nodes: Vec<Node>,
    cfg: PoaConfig,
    n_reads: usize,
    read_lens: Vec<usize>,
    /// Per-read provenance: the ordered node indices each read occupies
    /// (Match/Insert nodes, in read order; Deletes skip). Read 0 is the seed.
    /// Node indices are stable (nodes are only appended), so a path recorded
    /// when a read was added stays valid for the life of the graph. This drives
    /// per-partition rebuild membership and the coverage invariant.
    read_paths: Vec<Vec<usize>>,
    /// Per-read *matched* adjacencies: edge (a, b) is recorded only when the
    /// read reached `b` directly after `a` with NO intervening Delete — i.e. a
    /// genuine match/insert adjacency, excluding delete-bypass "resume" edges.
    /// This is the poa2 equivalent of legacy `edge_reads` (matched-only axis):
    /// bubble detection and arm-membership phasing use these, NOT the unified
    /// `Edge.weight` (which folds in delete-bypass traffic and would make a
    /// length variant's short-allele bypass look like a competing arm).
    read_matched_edges: Vec<Vec<(usize, usize)>>,
}

impl Poa {
    pub fn new(seed: &[u8], cfg: PoaConfig) -> Self {
        let mut nodes = Vec::with_capacity(seed.len());
        for (i, &b) in seed.iter().enumerate() {
            let inc = if i == 0 { vec![] } else { vec![i - 1] };
            let out = if i + 1 < seed.len() {
                vec![Edge {
                    to: i + 1,
                    weight: 1,
                }]
            } else {
                vec![]
            };
            nodes.push(Node {
                base: b,
                out,
                inc,
                aligned: vec![],
                cov: 1,
                del: 0,
            });
        }
        let seed_path: Vec<usize> = (0..seed.len()).collect();
        let seed_medges: Vec<(usize, usize)> = (0..seed.len().saturating_sub(1))
            .map(|i| (i, i + 1))
            .collect();
        Poa {
            nodes,
            cfg,
            n_reads: 1,
            read_lens: vec![seed.len()],
            read_paths: vec![seed_path],
            read_matched_edges: vec![seed_medges],
        }
    }

    /// Kahn's topological order over the current DAG. Returns (topo, rank_of).
    fn topo_order(&self) -> (Vec<usize>, Vec<usize>) {
        let n = self.nodes.len();
        let mut indeg = vec![0usize; n];
        for nd in &self.nodes {
            for e in &nd.out {
                indeg[e.to] += 1;
            }
        }
        let mut stack: Vec<usize> = (0..n).filter(|&i| indeg[i] == 0).collect();
        stack.sort_unstable(); // deterministic
        let mut topo = Vec::with_capacity(n);
        let mut rank_of = vec![0usize; n];
        while let Some(u) = stack.pop() {
            rank_of[u] = topo.len();
            topo.push(u);
            for e in &self.nodes[u].out {
                indeg[e.to] -= 1;
                if indeg[e.to] == 0 {
                    stack.push(e.to);
                }
            }
        }
        debug_assert_eq!(topo.len(), n, "graph must be a DAG");
        (topo, rank_of)
    }

    /// Full partial-order affine-gap GLOBAL alignment of `read` to the graph.
    /// Standard Gotoh with three states, generalized to a DAG: node `t`'s M/D
    /// relax over all predecessors. Returns ops in forward read order.
    fn align(&self, read: &[u8], topo: &[usize], rank_of: &[usize]) -> Vec<Op> {
        const VSTART: u32 = u32::MAX; // pred-rank sentinel: virtual start
        let n = self.nodes.len();
        let l = read.len();
        let (mm, mis, go, ge) = (
            self.cfg.match_score,
            self.cfg.mismatch_score,
            self.cfg.gap_open,
            self.cfg.gap_extend,
        );
        let semi = self.cfg.alignment_mode == AlignmentMode::SemiGlobal;

        // ---- Adaptive static-diagonal-union band (abPOA-style anti-fold) ----
        // Per row we score only query columns [lo, hi] = union of the
        // predecessor-following adaptive center with the static graph-geometry
        // diagonal, padded by `half`. The static diagonal is fixed by graph
        // geometry so it is ALWAYS in-band no matter how far the adaptive center
        // drifts -- this is what stops a repeat fold (abPOA GET_AD_DP_*).
        // Storage is BANDED: each row keeps only its [lo, hi] window in flat
        // arrays with per-row offsets, so memory is O(nodes × band) not
        // O(nodes × readlen).
        let banded0 = self.cfg.band_width > 0 || self.cfg.adaptive_band;
        let half = {
            let adaptive =
                self.cfg.adaptive_band_b + (self.cfg.adaptive_band_f * l as f32).ceil() as usize;
            adaptive.max(self.cfg.band_width).max(4)
        };
        let mut remain = vec![0usize; n]; // heaviest-out-edge node-distance to sink
        for &t in topo.iter().rev() {
            let mut best_w = -1i32;
            let mut best_c: Option<usize> = None;
            for e in &self.nodes[t].out {
                if e.weight > best_w {
                    best_w = e.weight;
                    best_c = Some(e.to);
                }
            }
            remain[t] = best_c.map_or(0, |c| remain[c] + 1);
        }
        // Band only spanning reads: `l - remain` assumes the read covers the
        // graph (false for a short partial read). Partials are short -> unbanded.
        let graph_span = remain.iter().max().copied().unwrap_or(0) + 1;
        let banded = banded0 && (l as f64) >= 0.8 * (graph_span as f64);

        // Banded per-row storage (indexed by rank). off[r]..off[r]+width is the
        // row's window; lo_[r]/hi_[r] its query-column bounds. j=0 column (D only,
        // leading deletions) is stored separately in d0/d0bk.
        let mut lo_ = vec![1usize; n];
        let mut hi_ = vec![0usize; n];
        let mut off = vec![0usize; n];
        let cap = if banded {
            n * (2 * half + 8).min(l + 1)
        } else {
            n * (l + 1)
        };
        let mut mflat: Vec<i32> = Vec::with_capacity(cap);
        let mut iflat: Vec<i32> = Vec::with_capacity(cap);
        let mut dflat: Vec<i32> = Vec::with_capacity(cap);
        let mut mbk: Vec<(u8, u32)> = Vec::with_capacity(cap);
        let mut dbk: Vec<(u8, u32)> = Vec::with_capacity(cap);
        let mut d0 = vec![NEG_INF; n];
        let mut d0bk = vec![(0u8, VSTART); n];
        let mut best_j = vec![1usize; n]; // per-rank winning column (adaptive center)

        for &t in topo.iter() {
            let r = rank_of[t];
            let base = self.nodes[t].base;
            let is_source = self.nodes[t].inc.is_empty();
            let (lo, hi) = if banded {
                let adaptive = if is_source {
                    1
                } else {
                    self.nodes[t]
                        .inc
                        .iter()
                        .map(|&p| best_j[rank_of[p]] + 1)
                        .max()
                        .unwrap_or(1)
                };
                let stat = l.saturating_sub(remain[t]).clamp(1, l.max(1));
                (
                    adaptive.min(stat).saturating_sub(half).max(1),
                    (adaptive.max(stat) + half).min(l).max(1),
                )
            } else {
                (1, l.max(1))
            };
            let hi = hi.max(lo);
            let width = hi + 1 - lo;
            lo_[r] = lo;
            hi_[r] = hi;
            off[r] = mflat.len();

            // ---- j = 0 column: D0 only (leading deletion chain) ----
            {
                let mut best = NEG_INF;
                let mut bk = (0u8, VSTART);
                if is_source {
                    // leading deletion from the virtual start
                    if go + ge > best {
                        best = go + ge;
                        bk = (0, VSTART);
                    }
                }
                for &p in &self.nodes[t].inc {
                    let pr = rank_of[p] as u32;
                    let de = d0[rank_of[p]];
                    if de != NEG_INF && de + ge > best {
                        best = de + ge;
                        bk = (2, pr);
                    }
                }
                d0[r] = best;
                d0bk[r] = bk;
            }

            let mut mrow = vec![NEG_INF; width];
            let mut irow = vec![NEG_INF; width];
            let mut drow = vec![NEG_INF; width];
            let mut mbrow = vec![(0u8, VSTART); width];
            let mut dbrow = vec![(0u8, VSTART); width];
            let mut row_best = NEG_INF;
            let mut row_best_j = best_j[r];

            for j in lo..=hi {
                let k = j - lo;
                // ---- M[t][j] ----
                let sc = if base == read[j - 1] { mm } else { mis };
                let mut best = NEG_INF;
                let mut bk = (0u8, VSTART);
                if (is_source || semi) && j == 1 && sc > best {
                    // free start: read[0] begins the alignment here
                    best = sc;
                    bk = (0, VSTART);
                }
                for &p in &self.nodes[t].inc {
                    let pr = rank_of[p];
                    let vm = bget(&mflat, &off, &lo_, &hi_, pr, j - 1);
                    let vi = bget(&iflat, &off, &lo_, &hi_, pr, j - 1);
                    let vd = bgetd(&dflat, &d0, &off, &lo_, &hi_, pr, j - 1);
                    for (st, v) in [(0u8, vm), (1u8, vi), (2u8, vd)] {
                        if v != NEG_INF && v + sc > best {
                            best = v + sc;
                            bk = (st, pr as u32);
                        }
                    }
                }
                mrow[k] = best;
                mbrow[k] = bk;
                // ---- I[t][j] ---- (read base inserted; same row, j-1) ----
                {
                    let (mo, io) = if j > lo {
                        (mrow[k - 1], irow[k - 1])
                    } else {
                        (NEG_INF, NEG_INF)
                    };
                    let mut best = NEG_INF;
                    if mo != NEG_INF {
                        best = mo + go + ge;
                    }
                    if io != NEG_INF && io + ge > best {
                        best = io + ge;
                    }
                    irow[k] = best;
                }
                // ---- D[t][j] ---- (node t skipped; predecessors at same j) ----
                {
                    let mut best = NEG_INF;
                    let mut bk = (0u8, VSTART);
                    for &p in &self.nodes[t].inc {
                        let pr = rank_of[p];
                        let mo = bget(&mflat, &off, &lo_, &hi_, pr, j);
                        let io = bget(&iflat, &off, &lo_, &hi_, pr, j);
                        let de = bgetd(&dflat, &d0, &off, &lo_, &hi_, pr, j);
                        if mo != NEG_INF && mo + go + ge > best {
                            best = mo + go + ge;
                            bk = (0, pr as u32);
                        }
                        if io != NEG_INF && io + go + ge > best {
                            best = io + go + ge;
                            bk = (1, pr as u32);
                        }
                        if de != NEG_INF && de + ge > best {
                            best = de + ge;
                            bk = (2, pr as u32);
                        }
                    }
                    drow[k] = best;
                    dbrow[k] = bk;
                }
                let cell = mrow[k].max(irow[k]).max(drow[k]);
                if cell != NEG_INF && cell > row_best {
                    row_best = cell;
                    row_best_j = j;
                }
            }
            best_j[r] = row_best_j;
            mflat.extend_from_slice(&mrow);
            iflat.extend_from_slice(&irow);
            dflat.extend_from_slice(&drow);
            mbk.extend_from_slice(&mbrow);
            dbk.extend_from_slice(&dbrow);
        }

        // Terminal (read fully consumed at j=l). Global: end at a sink; semi:
        // end at ANY node (free trailing graph gap).
        let mut best = NEG_INF;
        let mut cur: (u8, u32) = (0, VSTART);
        for &t in topo.iter() {
            if !semi && !self.nodes[t].out.is_empty() {
                continue;
            }
            let r = rank_of[t];
            let vm = bget(&mflat, &off, &lo_, &hi_, r, l);
            let vi = bget(&iflat, &off, &lo_, &hi_, r, l);
            let vd = bgetd(&dflat, &d0, &off, &lo_, &hi_, r, l);
            for (st, v) in [(0u8, vm), (1u8, vi), (2u8, vd)] {
                if v != NEG_INF && v > best {
                    best = v;
                    cur = (st, r as u32);
                }
            }
        }

        // Traceback.
        let mut ops: Vec<Op> = Vec::new();
        let (mut state, mut rr) = cur;
        let mut j = l;
        while rr != VSTART {
            let r = rr as usize;
            let t = topo[r];
            match state {
                0 => {
                    ops.push(Op::Match(t));
                    let (pst, pr) = mbk[off[r] + (j - lo_[r])];
                    j -= 1;
                    state = pst;
                    rr = pr;
                }
                1 => {
                    ops.push(Op::Ins(read[j - 1]));
                    j -= 1;
                    if j == 0 {
                        rr = VSTART;
                    } else {
                        let o = bget(&mflat, &off, &lo_, &hi_, r, j);
                        let e = bget(&iflat, &off, &lo_, &hi_, r, j);
                        state = if o >= e { 0 } else { 1 };
                    }
                }
                _ => {
                    ops.push(Op::Del(t));
                    let (pst, pr) = if j == 0 {
                        d0bk[r]
                    } else {
                        dbk[off[r] + (j - lo_[r])]
                    };
                    state = pst;
                    rr = pr;
                }
            }
            if rr == VSTART {
                while j > 0 {
                    ops.push(Op::Ins(read[j - 1]));
                    j -= 1;
                }
                break;
            }
        }
        ops.reverse();
        ops
    }

    fn add_edge(&mut self, from: usize, to: usize) {
        for e in &mut self.nodes[from].out {
            if e.to == to {
                e.weight += 1;
                return;
            }
        }
        self.nodes[from].out.push(Edge { to, weight: 1 });
        self.nodes[to].inc.push(from);
    }

    fn new_node(&mut self, base: u8) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(Node {
            base,
            out: vec![],
            inc: vec![],
            aligned: vec![],
            cov: 0,
            del: 0,
        });
        idx
    }

    /// SPOA/abPOA-style node fusion: the read base `rb` aligned (by the DP) to
    /// column-node `t`. If `t` (or one of its aligned siblings) already carries
    /// `rb`, reuse that node; otherwise create a new substitution node and
    /// mutually cross-link it into `t`'s alignment column. This keeps all
    /// alternative bases at a position as siblings, so later reads reuse them
    /// (correct allele counting) instead of fragmenting into fresh nodes.
    fn aligned_or_new(&mut self, t: usize, rb: u8) -> usize {
        if self.nodes[t].base == rb {
            return t;
        }
        for &a in &self.nodes[t].aligned {
            if self.nodes[a].base == rb {
                return a;
            }
        }
        let nn = self.new_node(rb);
        let mut column = self.nodes[t].aligned.clone();
        column.push(t);
        for &c in &column {
            self.nodes[c].aligned.push(nn);
            self.nodes[nn].aligned.push(c);
        }
        nn
    }

    /// Integrate a read's alignment into the graph and return `(path, medges)`:
    /// the ordered node indices it occupies (Match/Insert nodes, Deletes
    /// skipped) and its *matched* adjacencies (edges reached with no intervening
    /// Delete — the delete-bypass "resume" edge is excluded). The caller stores
    /// these in `read_paths` / `read_matched_edges`.
    fn integrate(&mut self, read: &[u8], ops: &[Op]) -> (Vec<usize>, Vec<(usize, usize)>) {
        // `prev` = last graph node the read is currently attached to.
        // Walk read positions: Match and Ins each consume one read base; Del none.
        let mut prev: Option<usize> = None;
        let mut j = 0usize;
        let mut path = Vec::with_capacity(read.len());
        let mut medges = Vec::with_capacity(read.len());
        // Whether a Delete occurred since the last occupied node: if so the next
        // Match/Ins reconnects via a bypass edge, which is NOT a matched
        // adjacency (mirrors legacy pure-bypass resume not touching edge_reads).
        let mut bypassed = false;
        for op in ops {
            match *op {
                Op::Match(t) => {
                    let rb = read[j];
                    j += 1;
                    let node = self.aligned_or_new(t, rb); // reuse on match, sub-node on mismatch
                    self.nodes[node].cov += 1;
                    if let Some(p) = prev {
                        self.add_edge(p, node);
                        if !bypassed {
                            medges.push((p, node));
                        }
                    }
                    prev = Some(node);
                    path.push(node);
                    bypassed = false;
                }
                Op::Ins(_) => {
                    let rb = read[j];
                    j += 1;
                    let nn = self.new_node(rb);
                    self.nodes[nn].cov += 1;
                    if let Some(p) = prev {
                        self.add_edge(p, nn);
                        if !bypassed {
                            medges.push((p, nn));
                        }
                    }
                    prev = Some(nn);
                    path.push(nn);
                    bypassed = false;
                }
                Op::Del(t) => {
                    // node skipped by this read: record the delete (used by the
                    // analysis layer's Match-vs-Delete column entropy) and mark
                    // that the next reconnection is a bypass, not a match.
                    self.nodes[t].del += 1;
                    bypassed = true;
                }
            }
        }
        (path, medges)
    }

    pub fn add_read(&mut self, read: &[u8]) {
        if read.is_empty() {
            return;
        }
        let (topo, rank_of) = self.topo_order();
        let ops = self.align(read, &topo, &rank_of);
        let (path, medges) = self.integrate(read, &ops);
        self.read_paths.push(path);
        self.read_matched_edges.push(medges);
        self.read_lens.push(read.len());
        self.n_reads += 1;
    }

    /// abPOA-style heaviest bundling: reverse pass computes, per node, the
    /// heaviest out-edge (tie-broken by downstream cumulative weight); walk the
    /// chosen chain from a source.
    /// abPOA-style heaviest bundling: the consensus node path (indices).
    fn heaviest_path_nodes(&self) -> Vec<usize> {
        let (topo, _) = self.topo_order();
        let n = self.nodes.len();
        if n == 0 {
            return vec![];
        }
        let mut score = vec![0i64; n];
        let mut nxt = vec![usize::MAX; n];
        for &t in topo.iter().rev() {
            let mut best_w = -1i32;
            let mut best_s = i64::MIN;
            let mut best_c = usize::MAX;
            for e in &self.nodes[t].out {
                let s = score[e.to];
                if e.weight > best_w || (e.weight == best_w && s > best_s) {
                    best_w = e.weight;
                    best_s = s;
                    best_c = e.to;
                }
            }
            if best_c != usize::MAX {
                score[t] = best_w as i64 + best_s;
                nxt[t] = best_c;
            }
        }
        let mut start = topo[0];
        let mut best = i64::MIN;
        for &t in &topo {
            if self.nodes[t].inc.is_empty() && score[t] > best {
                best = score[t];
                start = t;
            }
        }
        let mut path = Vec::new();
        let mut cur = start;
        while cur != usize::MAX {
            path.push(cur);
            cur = nxt[cur];
        }
        // Boundary trim: drop leading/trailing consensus nodes below the
        // coverage floor. This removes low-coverage flanks and trailing/leading
        // repeat-unit extensions supported by only a minority of reads (e.g. a
        // few longer reads over-extending a homogeneous repeat). abPOA/SPOA do
        // the equivalent; unlike the legacy engine we do NOT also run an
        // interior filter -- the band prevents the fold that motivated it.
        // Interior low-coverage (a genuine spanning-read gap) is preserved.
        let floor = coverage_threshold(self.n_reads, self.cfg.min_coverage_fraction);
        let s = path.iter().position(|&nd| self.nodes[nd].cov >= floor);
        let e = path.iter().rposition(|&nd| self.nodes[nd].cov >= floor);
        match (s, e) {
            (Some(s), Some(e)) if s <= e => path[s..=e].to_vec(),
            _ => path,
        }
    }

    pub fn consensus(&self) -> Vec<u8> {
        self.heaviest_path_nodes()
            .iter()
            .map(|&nd| self.nodes[nd].base)
            .collect()
    }

    /// Materialise an arm from `start` until the first reconvergence (2+ in-edges)
    /// or fork or dead end, capped at `cap` bases. Returns (sequence, node_count).
    fn materialize_arm(&self, start: usize, cap: usize) -> (Vec<u8>, usize) {
        let mut seq = Vec::new();
        let mut cur = start;
        let mut count = 0usize;
        loop {
            seq.push(self.nodes[cur].base);
            count += 1;
            if seq.len() >= cap {
                break;
            }
            match self.nodes[cur].out.as_slice() {
                [e] => {
                    let nxt = e.to;
                    if self.nodes[nxt].inc.len() > 1 {
                        break; // reconvergence
                    }
                    cur = nxt;
                }
                _ => break, // fork or dead end
            }
        }
        (seq, count)
    }

    /// GraphStats over the whole graph (O(V+E)).
    fn compute_stats(&self) -> crate::types::GraphStats {
        use crate::types::GraphStats;
        let n = self.nodes.len();
        let edge_count: usize = self.nodes.iter().map(|nd| nd.out.len()).sum();
        // coverage mean/variance + single-support fraction
        let mut sum = 0.0f64;
        let mut single = 0usize;
        for nd in &self.nodes {
            sum += nd.cov as f64;
            if nd.cov == 1 {
                single += 1;
            }
        }
        let mean = if n > 0 { sum / n as f64 } else { 0.0 };
        let mut var = 0.0f64;
        for nd in &self.nodes {
            let d = nd.cov as f64 - mean;
            var += d * d;
        }
        var = if n > 0 { var / n as f64 } else { 0.0 };
        // edge-weight Gini
        let mut ws: Vec<i32> = self
            .nodes
            .iter()
            .flat_map(|nd| nd.out.iter().map(|e| e.weight))
            .collect();
        ws.sort_unstable();
        let gini = if ws.is_empty() {
            0.0
        } else {
            let tot: i64 = ws.iter().map(|&w| w as i64).sum();
            if tot == 0 {
                0.0
            } else {
                let mut cum = 0i64;
                let mut area = 0i64;
                for &w in &ws {
                    cum += w as i64;
                    area += cum;
                }
                let nlen = ws.len() as f64;
                1.0 - (2.0 * area as f64) / (nlen * tot as f64) - 1.0 / nlen
            }
        };
        // bubbles: fork nodes with 2+ out-edges each >= allele floor
        let floor = ((self.n_reads as f64) * self.cfg.min_allele_freq).floor() as i32;
        let floor = floor.max(1);
        let mut bubble_count = 0usize;
        let mut max_bubble_depth = 0usize;
        let mut longest_span = 0usize;
        for nd in &self.nodes {
            let sig: Vec<i32> = nd
                .out
                .iter()
                .map(|e| e.weight)
                .filter(|&w| w >= floor)
                .collect();
            if sig.len() >= 2 {
                bubble_count += 1;
                let mut s = sig.clone();
                s.sort_unstable_by(|a, b| b.cmp(a));
                max_bubble_depth = max_bubble_depth.max(s[1] as usize); // minority arm
                for e in &nd.out {
                    if e.weight >= floor {
                        let (_, cnt) = self.materialize_arm(e.to, 4096);
                        longest_span = longest_span.max(cnt);
                    }
                }
            }
        }
        // mean column entropy over aligned columns
        let mut visited = vec![false; n];
        let mut ent_sum = 0.0f64;
        let mut ncol = 0usize;
        for i in 0..n {
            if visited[i] {
                continue;
            }
            let mut col = vec![i];
            col.extend(self.nodes[i].aligned.iter().copied());
            // Column distribution = each aligned base's coverage PLUS a "delete"
            // alternative (reads that skip the whole column) = the max
            // delete_count in the column. This makes mean_column_entropy measure
            // Match-vs-Delete (indel) disagreement as well as base substitution.
            let mut counts: Vec<f64> = col.iter().map(|&c| self.nodes[c].cov as f64).collect();
            let del = col.iter().map(|&c| self.nodes[c].del).max().unwrap_or(0);
            if del > 0 {
                counts.push(del as f64);
            }
            for &c in &col {
                visited[c] = true;
            }
            let tot: f64 = counts.iter().sum();
            if tot > 0.0 {
                let mut h = 0.0f64;
                for &cnt in &counts {
                    let p = cnt / tot;
                    if p > 0.0 {
                        h -= p * p.log2();
                    }
                }
                ent_sum += h;
                ncol += 1;
            }
        }
        let mean_column_entropy = if ncol > 0 { ent_sum / ncol as f64 } else { 0.0 };
        let mut lens = self.read_lens.clone();
        lens.sort_unstable();
        let median_input_read_len = if lens.is_empty() {
            0
        } else {
            lens[lens.len() / 2]
        };
        GraphStats {
            node_count: n,
            edge_count,
            bubble_count,
            max_bubble_depth,
            coverage_mean: mean,
            coverage_variance: var,
            edge_weight_gini: gini,
            single_support_fraction: if n > 0 { single as f64 / n as f64 } else { 0.0 },
            mean_column_entropy,
            longest_bubble_span: longest_span,
            median_input_read_len,
        }
    }

    /// Full consensus output (rich `Consensus`) for the public API.
    pub fn consensus_full(&self) -> crate::types::Consensus {
        use crate::types::{BubbleSite, Consensus, CoverageGap, GapKind};
        let path = self.heaviest_path_nodes();
        let sequence: Vec<u8> = path.iter().map(|&nd| self.nodes[nd].base).collect();
        let coverage: Vec<u32> = path.iter().map(|&nd| self.nodes[nd].cov).collect();
        let mut path_weights = Vec::with_capacity(path.len());
        for (i, &nd) in path.iter().enumerate() {
            if i == 0 {
                path_weights.push(self.nodes[nd].cov as i32);
            } else {
                let prev = path[i - 1];
                let w = self.nodes[prev]
                    .out
                    .iter()
                    .find(|e| e.to == nd)
                    .map(|e| e.weight)
                    .unwrap_or(0);
                path_weights.push(w);
            }
        }
        // interior coverage gaps: maximal runs at seed-only depth (cov <= 1)
        // flanked by higher coverage on both sides.
        let mut gaps = Vec::new();
        let cl = coverage.len();
        let mut i = 0;
        while i < cl {
            if coverage[i] <= 1 {
                let s = i;
                while i < cl && coverage[i] <= 1 {
                    i += 1;
                }
                if s > 0 && i < cl {
                    gaps.push(CoverageGap {
                        start: s,
                        end: i,
                        kind: GapKind::Spanning,
                    });
                }
            } else {
                i += 1;
            }
        }
        // bubble sites on the consensus path
        let floor = (((self.n_reads as f64) * self.cfg.min_allele_freq).floor() as i32).max(1);
        let mut bubble_sites = Vec::new();
        for (pos, &nd) in path.iter().enumerate() {
            let sig: Vec<&Edge> = self.nodes[nd]
                .out
                .iter()
                .filter(|e| e.weight >= floor)
                .collect();
            if sig.len() >= 2 {
                let mut arm_read_counts = Vec::new();
                let mut arm_sequences = Vec::new();
                let mut is_structural = false;
                for e in &sig {
                    arm_read_counts.push(e.weight as u32);
                    let (seq, cnt) = self.materialize_arm(e.to, 256);
                    if cnt >= self.cfg.phasing_bubble_min_span {
                        is_structural = true;
                    }
                    arm_sequences.push(if seq.len() >= 256 { Vec::new() } else { seq });
                }
                bubble_sites.push(BubbleSite {
                    consensus_pos: pos,
                    arm_read_counts,
                    arm_sequences,
                    is_structural,
                });
            }
        }
        Consensus {
            sequence,
            coverage,
            path_weights,
            n_reads: self.n_reads,
            graph_stats: self.compute_stats(),
            gaps,
            bubble_sites,
            read_indices: vec![],
        }
    }

    /// Per-edge read provenance on the *matched* axis (the legacy `edge_reads`
    /// equivalent): edge (a, b) is confirmed by every read that reached `b`
    /// directly after `a` with no intervening Delete. Built from
    /// `read_matched_edges`, so delete-bypass resume edges are excluded — a
    /// length variant's short-allele bypass does NOT register as arm membership
    /// (matching legacy, which compensates via read-length routing).
    fn edge_reads(&self) -> HashMap<(usize, usize), Vec<u32>> {
        let mut m: HashMap<(usize, usize), Vec<u32>> = HashMap::new();
        for (r, medges) in self.read_matched_edges.iter().enumerate() {
            for &e in medges {
                m.entry(e).or_default().push(r as u32);
            }
        }
        for v in m.values_mut() {
            v.sort_unstable();
            v.dedup();
        }
        m
    }

    /// Build the read × het-site call matrix for linkage phasing. A het site is
    /// a fork (a node with ≥2 out-edges each clearing `floor` reads); node fusion
    /// makes substitutions forks too, so this uniformly captures SNP, indel, and
    /// length divergences. A read's call at a site is the branch (successor) its
    /// matched path takes, or MISSING if it doesn't traverse an out-edge of the
    /// fork or takes a below-floor minor branch. Uses `read_matched_edges`
    /// (matched-only, excludes delete-bypass resume) so a length variant's short
    /// allele registers as a genuine branch, consistent with the phasing signal.
    pub fn phasing_matrix(&self, floor: u32) -> crate::phasing::PhasingMatrix {
        use crate::phasing::{MISSING, PhasingMatrix, Site};
        // Per-read matched adjacency map: node -> successor taken by this read.
        let mut read_next: Vec<std::collections::HashMap<usize, usize>> =
            vec![std::collections::HashMap::new(); self.n_reads];
        for (r, medges) in self.read_matched_edges.iter().enumerate() {
            for &(a, b) in medges {
                read_next[r].insert(a, b);
            }
        }
        // Matched-only out-edge weights: (from,to) -> reads.
        let mut mw: std::collections::HashMap<(usize, usize), u32> =
            std::collections::HashMap::new();
        for medges in &self.read_matched_edges {
            for &e in medges {
                *mw.entry(e).or_default() += 1;
            }
        }
        // Fork sites: nodes with ≥2 matched out-edges each ≥ floor.
        let mut sites: Vec<Site> = Vec::new();
        for node in 0..self.nodes.len() {
            let mut branches: Vec<(usize, u32)> = self.nodes[node]
                .out
                .iter()
                .filter_map(|e| {
                    let w = mw.get(&(node, e.to)).copied().unwrap_or(0);
                    if w >= floor { Some((e.to, w)) } else { None }
                })
                .collect();
            if branches.len() < 2 {
                continue;
            }
            branches.sort_unstable_by_key(|&(to, _)| to); // deterministic
            sites.push(Site {
                node,
                branches: branches.iter().map(|&(to, _)| to).collect(),
                branch_support: branches.iter().map(|&(_, w)| w).collect(),
            });
        }
        // Calls: for each read/site, which branch (successor) did it take.
        let mut calls = vec![vec![MISSING; sites.len()]; self.n_reads];
        for (s, site) in sites.iter().enumerate() {
            for r in 0..self.n_reads {
                if let Some(&nxt) = read_next[r].get(&site.node) {
                    if let Some(bi) = site.branches.iter().position(|&b| b == nxt) {
                        calls[r][s] = bi as i16;
                    }
                }
            }
        }
        PhasingMatrix {
            n_reads: self.n_reads,
            sites,
            calls,
        }
    }

    /// Matched-only adjacency view for bubble detection: `(out, inc)` where
    /// `out[a]` is `[(b, matched_weight)]` and `inc[b]` is `[a, …]`, both
    /// counting only matched adjacencies (delete-bypass excluded). This is the
    /// legacy `Node.out_edges`/`ew.matched` + `in_edges` semantics, rebuilt so
    /// the ported bubble finders never see poa2's bypass-inclusive `Edge.weight`.
    fn matched_view(&self) -> (MatchedOut, MatchedInc) {
        let n = self.nodes.len();
        let mut w: HashMap<(usize, usize), i32> = HashMap::new();
        for medges in &self.read_matched_edges {
            for &e in medges {
                *w.entry(e).or_default() += 1;
            }
        }
        let mut out = vec![Vec::new(); n];
        let mut inc = vec![Vec::new(); n];
        // Deterministic ordering: emit each node's out-edges in ascending `to`.
        let mut keys: Vec<(usize, usize)> = w.keys().copied().collect();
        keys.sort_unstable();
        for (a, b) in keys {
            out[a].push((b, w[&(a, b)]));
            inc[b].push(a);
        }
        (out, inc)
    }
}

// ─── Multi-allele extraction (ported from graph.rs) ───────────────────────────
//
// The pipeline is faithful to the legacy `graph.rs` implementation (see the
// crate CLAUDE.md "Multi-Allele Extraction" section for the full rationale
// behind every threshold and discriminator). The only adaptations are:
//   * the node model: `out_edges: Vec<(usize, EdgeWeight)>` / `ew.matched`
//     becomes poa2's `out: Vec<Edge>` / `e.weight`; `in_edges` becomes `inc`.
//   * provenance: legacy `edge_reads` is derived from `read_paths` (see
//     `Poa::edge_reads`); read indices are internal (graph add order).
//   * read lengths: the clustering helpers take a `lens: &[usize]` slice
//     (internal-indexed) instead of `&[Vec<u8>]`.
//   * per-partition consensus is a fresh poa2 rebuild + `consensus_full` (the
//     legacy path likewise rebuilds a sub-graph per group), so there is no
//     restricted-heaviest-path machinery.

/// Matched-only adjacency view: `out[a] = [(b, matched_weight)]`.
type MatchedOut = Vec<Vec<(usize, i32)>>;
/// Matched-only reverse adjacency: `inc[b] = [a, …]`.
type MatchedInc = Vec<Vec<usize>>;

const MIN_LENGTH_GAP_BP: f64 = 8.0;
const LENGTH_SEPARATION_MADS: f64 = 3.0;
const MIN_SPREAD_FLOOR_BP: f64 = 3.0;

fn median_usize(vals: &[usize]) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    let mut v: Vec<f64> = vals.iter().map(|&x| x as f64).collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

fn mad_usize(vals: &[usize], median: f64) -> f64 {
    let dev: Vec<usize> = vals
        .iter()
        .map(|&x| (x as f64 - median).abs().round() as usize)
        .collect();
    median_usize(&dev)
}

/// True when two groups' read-length distributions are genuinely bimodally
/// separated (real length-differing alleles) rather than ordinary scatter
/// around one shared mode.
fn length_separated(lens_a: &[usize], lens_b: &[usize]) -> bool {
    if lens_a.is_empty() || lens_b.is_empty() {
        return false;
    }
    let med_a = median_usize(lens_a);
    let med_b = median_usize(lens_b);
    let gap = (med_a - med_b).abs();
    if gap < MIN_LENGTH_GAP_BP {
        return false;
    }
    let spread_a = mad_usize(lens_a, med_a).max(MIN_SPREAD_FLOOR_BP);
    let spread_b = mad_usize(lens_b, med_b).max(MIN_SPREAD_FLOOR_BP);
    gap >= LENGTH_SEPARATION_MADS * spread_a.max(spread_b)
}

/// Number of `node`'s in-edges (matched-view) whose own source-side weight
/// clears `sig_threshold` (a genuine reconvergence vs one real path plus noise).
fn real_in_edge_count(
    out: &MatchedOut,
    inc: &MatchedInc,
    node: usize,
    sig_threshold: i32,
) -> usize {
    inc[node]
        .iter()
        .filter(|&&p| {
            out[p]
                .iter()
                .find(|&&(to, _)| to == node)
                .is_some_and(|&(_, w)| w >= sig_threshold)
        })
        .count()
}

/// Noise-tolerant arm-span walk over the matched-view: a fork or reconvergence
/// ends the arm only when it is a *real* one (2+ edges each clearing
/// `sig_threshold`). See the legacy `materialize_arm_len_tolerant` for why safe.
fn materialize_arm_len_tolerant(
    out: &MatchedOut,
    inc: &MatchedInc,
    start: usize,
    max_depth: usize,
    sig_threshold: i32,
) -> usize {
    let mut len = 0;
    let mut cur = start;
    for _ in 0..max_depth {
        if len >= 1 && real_in_edge_count(out, inc, cur, sig_threshold) > 1 {
            break;
        }
        len += 1;
        match out[cur].as_slice() {
            [] => break,
            [(next, _)] => cur = *next,
            edges => {
                let real_arms = edges.iter().filter(|&&(_, w)| w >= sig_threshold).count();
                if real_arms > 1 {
                    break; // a genuine nested fork -- this arm ends here.
                }
                let &(next, _) = edges.iter().max_by_key(|&&(_, w)| w).unwrap();
                cur = next;
            }
        }
    }
    len
}

/// Bubbles where at least one arm has structural size (span ≥
/// `phasing_bubble_min_span`): genuine length variants / SVs, not SNPs.
fn find_structural_bubbles(
    out: &MatchedOut,
    inc: &MatchedInc,
    topo: &[usize],
    n_reads: usize,
    cfg: &PoaConfig,
) -> Vec<(usize, Vec<usize>)> {
    let threshold = ((n_reads as f64 * cfg.min_allele_freq).ceil() as i32).max(1);
    let max_check = cfg.phasing_bubble_min_span.saturating_add(1);

    topo.iter()
        .copied()
        .filter_map(|entry_node| {
            let arms: Vec<usize> = out[entry_node]
                .iter()
                .filter(|&&(_, w)| w >= threshold)
                .map(|&(to, _)| to)
                .collect();
            if arms.len() < 2 {
                return None;
            }
            let max_span = arms
                .iter()
                .map(|&start| {
                    if real_in_edge_count(out, inc, start, threshold) > 1 {
                        0
                    } else {
                        materialize_arm_len_tolerant(out, inc, start, max_check, threshold)
                    }
                })
                .max()
                .unwrap_or(0);
            if max_span >= cfg.phasing_bubble_min_span {
                Some((entry_node, arms))
            } else {
                None
            }
        })
        .collect()
}

/// Bubbles with 2+ arms each clearing the `min_allele_freq` vote threshold
/// (any span; the SNP-fallback finder) over the matched-view.
fn find_bubbles(
    out: &MatchedOut,
    topo: &[usize],
    n_reads: usize,
    min_allele_freq: f64,
) -> Vec<(usize, Vec<usize>)> {
    let threshold = ((n_reads as f64 * min_allele_freq).ceil() as i32).max(1);
    topo.iter()
        .copied()
        .filter_map(|node_idx| {
            let arms: Vec<usize> = out[node_idx]
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

/// Groups reads by arm-choice compatibility across all structural bubbles.
/// Incremental, most-informative-first clustering (NOT plain union-find); bridge
/// candidates and unassigned reads resolved by length with an 85% plausibility
/// gate. `lens` is internal-read-indexed.
fn phasing_groups(
    edge_reads: &HashMap<(usize, usize), Vec<u32>>,
    bubbles: &[(usize, Vec<usize>)],
    n_reads: usize,
    min_reads: usize,
    lens: &[usize],
) -> Vec<Vec<usize>> {
    let n_bubbles = bubbles.len();

    let mut sig: Vec<Vec<Option<usize>>> = vec![vec![None; n_bubbles]; n_reads];
    for (b, (entry, arm_starts)) in bubbles.iter().enumerate() {
        for (arm_idx, &arm_start) in arm_starts.iter().enumerate() {
            if let Some(rs) = edge_reads.get(&(*entry, arm_start)) {
                for &r in rs {
                    let r = r as usize;
                    if r < n_reads {
                        sig[r][b] = Some(arm_idx);
                    }
                }
            }
        }
    }

    let mut assigned: Vec<usize> = Vec::new();
    let mut unassigned: Vec<usize> = Vec::new();
    for (r, row) in sig.iter().enumerate() {
        if row.iter().any(|s| s.is_some()) {
            assigned.push(r);
        } else {
            unassigned.push(r);
        }
    }

    let mut order = assigned.clone();
    order.sort_by_key(|&r| {
        let known = sig[r].iter().filter(|s| s.is_some()).count();
        (std::cmp::Reverse(known), r)
    });

    let mut clusters: Vec<Vec<Option<usize>>> = Vec::new();
    let mut cluster_members: Vec<Vec<usize>> = Vec::new();
    let mut bridge_candidates: Vec<(usize, Vec<usize>)> = Vec::new();

    for r in order {
        let row = &sig[r];
        let compatible_clusters: Vec<usize> = clusters
            .iter()
            .enumerate()
            .filter(|(_, known)| {
                (0..n_bubbles).all(|b| match (row[b], known[b]) {
                    (Some(a), Some(bv)) => a == bv,
                    _ => true,
                })
            })
            .map(|(ci, _)| ci)
            .collect();

        match compatible_clusters.as_slice() {
            [] => {
                clusters.push(row.clone());
                cluster_members.push(vec![r]);
            }
            [only] => {
                cluster_members[*only].push(r);
                for b in 0..n_bubbles {
                    if clusters[*only][b].is_none() {
                        clusters[*only][b] = row[b];
                    }
                }
            }
            _ => bridge_candidates.push((r, compatible_clusters)),
        }
    }

    const PLAUSIBLE_LEN_FRACTION: f64 = 0.85;

    if !cluster_members.is_empty() {
        let snapshot_lens: Vec<f64> = cluster_members
            .iter()
            .map(|members| median_usize(&members.iter().map(|&r| lens[r]).collect::<Vec<_>>()))
            .collect();
        let snapshot_sizes: Vec<usize> = cluster_members.iter().map(|m| m.len()).collect();
        let min_full_len = (0..cluster_members.len())
            .filter(|&c| snapshot_sizes[c] >= min_reads)
            .map(|c| snapshot_lens[c])
            .fold(f64::INFINITY, f64::min);

        let nearest_cluster_among = |read_idx: usize, candidates: &[usize]| -> usize {
            let len = lens[read_idx] as f64;
            candidates
                .iter()
                .copied()
                .min_by(|&a, &b| {
                    let da = (snapshot_lens[a] - len).abs();
                    let db = (snapshot_lens[b] - len).abs();
                    da.partial_cmp(&db)
                        .unwrap()
                        .then_with(|| snapshot_sizes[b].cmp(&snapshot_sizes[a]))
                        .then_with(|| a.cmp(&b))
                })
                .unwrap()
        };
        let largest_overall = (0..cluster_members.len())
            .max_by_key(|&c| snapshot_sizes[c])
            .unwrap();
        let all_clusters: Vec<usize> = (0..cluster_members.len()).collect();
        let largest_compatible = |candidates: &[usize]| -> usize {
            *candidates
                .iter()
                .max_by_key(|&&c| snapshot_sizes[c])
                .unwrap()
        };

        let mut resolved: Vec<(usize, usize)> = Vec::new();
        for (r, compat) in &bridge_candidates {
            let len = lens[*r] as f64;
            let target = if min_full_len.is_finite() && len >= PLAUSIBLE_LEN_FRACTION * min_full_len
            {
                nearest_cluster_among(*r, &all_clusters)
            } else {
                largest_compatible(compat)
            };
            resolved.push((*r, target));
        }
        for &r in &unassigned {
            let len = lens[r] as f64;
            let target = if min_full_len.is_finite() && len >= PLAUSIBLE_LEN_FRACTION * min_full_len
            {
                nearest_cluster_among(r, &all_clusters)
            } else {
                largest_overall
            };
            resolved.push((r, target));
        }
        for (r, target) in resolved {
            cluster_members[target].push(r);
        }
    } else {
        cluster_members.push(Vec::new());
        cluster_members[0].extend(bridge_candidates.into_iter().map(|(r, _)| r));
        cluster_members[0].extend(unassigned);
    }

    let mut groups: Vec<Vec<usize>> = cluster_members;
    groups.sort_unstable_by_key(|g| std::cmp::Reverse(g.len()));
    if groups.is_empty() {
        groups.push(Vec::new());
    }
    let mut i = 1;
    while i < groups.len() {
        if groups[i].len() < min_reads {
            let g = groups.remove(i);
            groups[0].extend(g);
        } else {
            i += 1;
        }
    }
    groups.retain(|g| !g.is_empty());
    groups
}

/// Safety net over `phasing_groups`: a provisional cluster is trusted as its
/// own allele only if it is length-separated from AND has a clean distinguishing
/// bubble vs every confirmed group (both gated on 2+ structural bubbles).
/// Rejected candidates route per-member to the length-nearest confirmed group.
fn validate_and_merge_groups(
    groups: Vec<Vec<usize>>,
    lens: &[usize],
    min_reads: usize,
    bubbles: &[(usize, Vec<usize>)],
    edge_reads: &HashMap<(usize, usize), Vec<u32>>,
) -> Vec<Vec<usize>> {
    if groups.len() < 2 {
        return groups;
    }
    let mut groups = groups;
    groups.sort_unstable_by_key(|g| std::cmp::Reverse(g.len()));

    let n_bubbles = bubbles.len();
    let require_length_separation = n_bubbles >= 2;

    const CLEAR_MAJORITY: f64 = 0.60;
    const MIN_ARM_COV: usize = 2;

    let arm_of: Vec<HashMap<usize, usize>> = bubbles
        .iter()
        .map(|(entry, arm_starts)| {
            let mut m: HashMap<usize, usize> = HashMap::new();
            for (k, &a) in arm_starts.iter().enumerate() {
                if let Some(rs) = edge_reads.get(&(*entry, a)) {
                    for &r in rs {
                        m.entry(r as usize).or_insert(k);
                    }
                }
            }
            m
        })
        .collect();

    let group_bubble_stats = |grp: &[usize], b: usize| -> (Option<usize>, f64, usize) {
        let mut counts: HashMap<usize, usize> = HashMap::new();
        for &r in grp {
            if let Some(&k) = arm_of[b].get(&r) {
                *counts.entry(k).or_default() += 1;
            }
        }
        let cov: usize = counts.values().sum();
        match counts.into_iter().max_by_key(|&(_, c)| c) {
            Some((arm, c)) => (Some(arm), c as f64 / cov as f64, cov),
            None => (None, 0.0, 0),
        }
    };

    let clean_distinguishing = |ga: &[usize], gb: &[usize]| -> bool {
        (0..n_bubbles).any(|b| {
            let (ma, fa, ca) = group_bubble_stats(ga, b);
            let (mb, fb, cb) = group_bubble_stats(gb, b);
            ma.is_some()
                && mb.is_some()
                && ma != mb
                && fa >= CLEAR_MAJORITY
                && fb >= CLEAR_MAJORITY
                && ca >= MIN_ARM_COV
                && cb >= MIN_ARM_COV
        })
    };

    let mut confirmed: Vec<Vec<usize>> = vec![groups[0].clone()];
    let mut confirmed_lens: Vec<Vec<usize>> = vec![groups[0].iter().map(|&r| lens[r]).collect()];

    for cand in groups.into_iter().skip(1) {
        let cand_lens: Vec<usize> = cand.iter().map(|&r| lens[r]).collect();

        let significant = cand.len() >= min_reads;
        let separated_from_all = !require_length_separation
            || confirmed_lens
                .iter()
                .all(|existing| length_separated(&cand_lens, existing));

        if significant && separated_from_all {
            let structurally_distinct = !require_length_separation
                || confirmed.iter().all(|c| clean_distinguishing(&cand, c));
            if structurally_distinct {
                confirmed.push(cand);
                confirmed_lens.push(cand_lens);
            } else {
                let cand_med = median_usize(&cand_lens);
                let target = (0..confirmed.len())
                    .filter(|&c| !clean_distinguishing(&cand, &confirmed[c]))
                    .min_by(|&a, &b| {
                        let da = (median_usize(&confirmed_lens[a]) - cand_med).abs();
                        let db = (median_usize(&confirmed_lens[b]) - cand_med).abs();
                        da.partial_cmp(&db).unwrap()
                    })
                    .expect("structurally_distinct false => at least one indistinguishable group");
                confirmed_lens[target].extend(cand_lens);
                confirmed[target].extend(cand);
            }
        } else if confirmed.len() == 1 {
            confirmed_lens[0].extend(cand_lens);
            confirmed[0].extend(cand);
        } else {
            let snapshot_medians: Vec<f64> =
                confirmed_lens.iter().map(|l| median_usize(l)).collect();
            let mut per_read_targets: Vec<(usize, usize)> = Vec::with_capacity(cand.len());
            for (&r, &len) in cand.iter().zip(cand_lens.iter()) {
                let len = len as f64;
                let target = snapshot_medians
                    .iter()
                    .enumerate()
                    .min_by(|&(_, &a), &(_, &b)| {
                        let da = (a - len).abs();
                        let db = (b - len).abs();
                        da.partial_cmp(&db).unwrap()
                    })
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                per_read_targets.push((r, target));
            }
            for (r, target) in per_read_targets {
                confirmed_lens[target].push(lens[r]);
                confirmed[target].push(r);
            }
        }
    }

    confirmed
}

/// SNP-fallback partition: assign each read to the first of the top-2 arms
/// (by read count) it appears in; unassigned reads go to the largest group.
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
        if let Some(rs) = edge_reads.get(&(entry_node, arm_start)) {
            for &r in rs {
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

/// Build a poa2 graph from `reads`, seeding on the median-length (non-empty)
/// read. Returns `(graph, ext)` where `ext[internal] = external index into
/// reads`; `graph.n_reads == ext.len()` (empty reads are skipped consistently).
fn build_median(reads: &[&[u8]], cfg: &PoaConfig) -> (Poa, Vec<usize>) {
    let mut order: Vec<usize> = (0..reads.len()).filter(|&i| !reads[i].is_empty()).collect();
    order.sort_by_key(|&i| reads[i].len());
    let med = order[order.len() / 2];
    let mut g = Poa::new(reads[med], cfg.clone());
    let mut ext = vec![med];
    for (i, r) in reads.iter().enumerate() {
        if i != med && !r.is_empty() {
            g.add_read(r);
            ext.push(i);
        }
    }
    (g, ext)
}

/// Classic O(nm) edit distance (small consensus sequences, exact).
fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let (n, m) = (a.len(), b.len());
    let mut prev: Vec<usize> = (0..=m).collect();
    let mut cur = vec![0usize; m + 1];
    for i in 1..=n {
        cur[0] = i;
        for j in 1..=m {
            let cost = usize::from(a[i - 1] != b[j - 1]);
            cur[j] = (prev[j] + 1).min(cur[j - 1] + 1).min(prev[j - 1] + cost);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[m]
}

/// Linkage-based multi-allele consensus (design/multi_allele_phasing_design.md).
///
/// Phases reads by **linkage** — a genuine allele is a set of differences that
/// co-occur on the same reads across many het sites — rather than by any single
/// bubble's arm span (a folded-graph artifact). The dominant read bipartition is
/// only *trusted* when its two sides yield genuinely different consensus
/// sequences (the consensus-difference confirmation), so phase-shift noise in a
/// homogeneous repeat — which splits reads but yields the SAME consensus on both
/// sides — collapses back to one allele. Recurses for 3+ alleles. Each allele's
/// `read_indices` are the external (input) read indices.
pub fn linkage_consensus_multi(
    reads: &[&[u8]],
    config: &PoaConfig,
) -> Result<Vec<crate::types::Consensus>, crate::error::PoaError> {
    use crate::error::PoaError;
    if reads.len() < config.min_reads {
        return Err(PoaError::InsufficientDepth {
            got: reads.len(),
            min: config.min_reads,
        });
    }
    let all: Vec<usize> = (0..reads.len()).filter(|&i| !reads[i].is_empty()).collect();
    let mut out = Vec::new();
    linkage_split(reads, &all, config, 0, &mut out);
    Ok(out)
}

/// Max recursion depth (allele count is bounded by 2^depth; 3 covers up to ~8).
const LINKAGE_MAX_DEPTH: usize = 3;
/// A read bipartition is only trusted as a real allele split when it is this
/// clean (fraction of calls matching their side's per-site majority). Noise
/// sub-clusters of one allele split at lower consistency, so this gate stops the
/// recursion from over-splitting a single allele at higher error rates.
const LINKAGE_SPLIT_CONSISTENCY: f64 = 0.80;
/// Minimum consensus-length gap (bp) between two sides for a *length*-variant
/// split. Above per-read stutter (±1 repeat unit) for the common motifs, so an
/// allele's own stutter sub-clusters don't get split, while a genuine
/// repeat-count difference (many bp) does.
const LINKAGE_LEN_GAP_BP: usize = 8;

/// Build a consensus for `subset` of `reads`, tagged with the subset's external
/// indices (sorted).
fn subset_consensus(
    reads: &[&[u8]],
    subset: &[usize],
    config: &PoaConfig,
) -> crate::types::Consensus {
    let sub_reads: Vec<&[u8]> = subset.iter().map(|&i| reads[i]).collect();
    let (g, _) = build_median(&sub_reads, config);
    let mut c = g.consensus_full();
    let mut idx = subset.to_vec();
    idx.sort_unstable();
    c.read_indices = idx;
    c
}

/// Recursively split `subset` by the dominant linkage bipartition, keeping a
/// split only when the two sides' consensuses genuinely differ.
fn linkage_split(
    reads: &[&[u8]],
    subset: &[usize],
    config: &PoaConfig,
    depth: usize,
    out: &mut Vec<crate::types::Consensus>,
) {
    // Too shallow to split into two min_reads-sized alleles, or too deep.
    if subset.len() < 2 * config.min_reads || depth >= LINKAGE_MAX_DEPTH {
        out.push(subset_consensus(reads, subset, config));
        return;
    }
    let sub_reads: Vec<&[u8]> = subset.iter().map(|&i| reads[i]).collect();
    let (g, ext) = build_median(&sub_reads, config);
    let floor = ((g.n_reads as f64 * config.min_allele_freq).ceil() as u32)
        .max(config.min_reads as u32)
        .max(2);
    let bp = g.phasing_matrix(floor).dominant_bipartition();

    // Cheap pre-filter: need a real candidate split (≥2 informative sites, both
    // sides big enough). The consensus-difference test below is the real gate.
    let mut s0 = Vec::new();
    let mut s1 = Vec::new();
    for internal in 0..g.n_reads {
        let orig = subset[ext[internal]];
        if bp.assign[internal] == 0 {
            s0.push(orig);
        } else {
            s1.push(orig);
        }
    }
    let whole = |out: &mut Vec<crate::types::Consensus>| {
        let mut c = g.consensus_full();
        let mut idx = subset.to_vec();
        idx.sort_unstable();
        c.read_indices = idx;
        out.push(c);
    };
    if bp.n_informative_sites < 2 || s0.len() < config.min_reads || s1.len() < config.min_reads {
        whole(out);
        return;
    }

    // Confirmation: do the two sides yield genuinely different sequences? A
    // phase-shift-noise split yields the same consensus on both sides; a real
    // allele split differs at the variant positions. Two regimes:
    //  - a LARGE difference (length variant) is real regardless of how noisy the
    //    per-read bipartition is (fixes under-splitting at high error), while
    //  - a SMALL difference (a few substitutions) is only trusted when the
    //    bipartition is clean (high consistency), so error noise between two
    //    sub-clusters of ONE allele does not fabricate a split.
    let c0 = subset_consensus(reads, &s0, config);
    let c1 = subset_consensus(reads, &s1, config);
    let diff = edit_distance(&c0.sequence, &c1.sequence);
    // A split is a real allele boundary in one of two regimes:
    //  - LENGTH variant: the two sides' consensus lengths differ by a real gap.
    //    Consensus length is the modal allele length, so per-read stutter and
    //    base errors (which don't shift the mode) leave this ~0 within one
    //    allele — this is what stops the recursion splitting an allele into its
    //    stutter sub-clusters at high error.
    //  - SUBSTITUTION haplotype: same length, but a clean (high-consistency)
    //    bipartition with base differences. Gated on consistency so error noise
    //    between two sub-clusters of one allele can't fabricate a split.
    let len_gap = c0.sequence.len().abs_diff(c1.sequence.len());
    let real_split = len_gap >= LINKAGE_LEN_GAP_BP
        || (bp.consistency >= LINKAGE_SPLIT_CONSISTENCY && len_gap <= 4 && diff >= 3);
    if !real_split {
        whole(out);
        return;
    }
    // Real split — recurse each side (handles 3+ alleles).
    linkage_split(reads, &s0, config, depth + 1, out);
    linkage_split(reads, &s1, config, depth + 1, out);
}

/// Multi-allele consensus over `reads` (WORK IN PROGRESS — not the production
/// path; `crate::consensus_multi` routes through the legacy engine). Structural
/// -bubble phasing (cross-bubble clustering + validate/merge), an unbanded
/// rebuild when no structural bubble is found on a banded graph, then a
/// single-best SNP-bubble fallback; each partition is rebuilt into its own poa2
/// sub-graph for the per-allele consensus.
///
/// KNOWN LIMITATION (does not yet reach legacy real-data parity): poa2's
/// node-fusion folds periodic repeats, so a shorter allele's length difference
/// is smeared across many equal-length phase bubbles instead of surfacing as one
/// high-span structural bubble the way legacy's separate `bypass_edges`
/// representation exposed it. Empirically, no phasing signal tried matches legacy
/// on the real HiFi multi scenarios: structural-only under-splits the clean
/// synthetic 3-allele case and fails multi_cag20_50 / multi_skew_cag20_40 /
/// multi_gaa30_100; all-bubbles over-splits multi_gaa30_100 (best at 19/20);
/// length-refinement over-splits worse (real reads have within-allele length
/// variance). Reaching parity needs flanking-anchor pre-processing
/// (`extract_repeat_segment`), a separate effort. This code carries the ported
/// clustering logic (correct given the right bubble set) for that follow-up.
pub fn consensus_multi(
    reads: &[&[u8]],
    config: &PoaConfig,
) -> Result<Vec<crate::types::Consensus>, crate::error::PoaError> {
    use crate::error::PoaError;
    if reads.len() < config.min_reads {
        return Err(PoaError::InsufficientDepth {
            got: reads.len(),
            min: config.min_reads,
        });
    }
    consensus_multi_impl(reads, config, true)
}

fn consensus_multi_impl(
    reads: &[&[u8]],
    config: &PoaConfig,
    allow_unbanded_retry: bool,
) -> Result<Vec<crate::types::Consensus>, crate::error::PoaError> {
    use crate::error::PoaError;

    let (g, ext) = build_median(reads, config);
    let lens: Vec<usize> = ext.iter().map(|&e| reads[e].len()).collect();
    let (topo, _) = g.topo_order();
    let (mout, min) = g.matched_view();

    let structural = find_structural_bubbles(&mout, &min, &topo, g.n_reads, config);

    // No structural bubble on a banded graph: the band may have folded a real
    // length variant onto the spine. Rebuild once fully unbanded and retry.
    if structural.is_empty()
        && allow_unbanded_retry
        && (config.band_width > 0 || config.adaptive_band)
    {
        let mut cfg2 = config.clone();
        cfg2.band_width = 0;
        cfg2.adaptive_band = false;
        return consensus_multi_impl(reads, &cfg2, false);
    }

    let edge_reads = g.edge_reads();

    let groups: Vec<Vec<usize>> = if !structural.is_empty() {
        let g0 = phasing_groups(&edge_reads, &structural, g.n_reads, config.min_reads, &lens);
        validate_and_merge_groups(g0, &lens, config.min_reads, &structural, &edge_reads)
    } else {
        let bubbles = find_bubbles(&mout, &topo, g.n_reads, config.min_allele_freq);
        if bubbles.is_empty() {
            return Ok(vec![g.consensus_full()]);
        }
        // Choose the bubble whose weakest arm has the most read support.
        let (entry, arms) = bubbles
            .iter()
            .max_by_key(|(entry, arms)| {
                arms.iter()
                    .map(|&a| edge_reads.get(&(*entry, a)).map_or(0, |v| v.len()))
                    .min()
                    .unwrap_or(0)
            })
            .unwrap();
        partition_reads_by_bubble(&edge_reads, *entry, arms, g.n_reads)
    };

    if groups.len() < 2 {
        return Ok(vec![g.consensus_full()]);
    }

    let mut out = Vec::with_capacity(groups.len());
    for group in groups {
        if group.len() < config.min_reads {
            return Err(PoaError::InsufficientDepth {
                got: group.len(),
                min: config.min_reads,
            });
        }
        let mut ext_idx: Vec<usize> = group.iter().map(|&i| ext[i]).collect();
        let sub_reads: Vec<&[u8]> = ext_idx.iter().map(|&e| reads[e]).collect();
        let mut sub_cfg = config.clone();
        sub_cfg.multi_allele = false;
        let (sub, _) = build_median(&sub_reads, &sub_cfg);
        let mut c = sub.consensus_full();
        ext_idx.sort_unstable();
        c.read_indices = ext_idx;
        out.push(c);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg() -> PoaConfig {
        PoaConfig::default()
    }

    #[test]
    fn identical_reads_roundtrip() {
        let seed = b"ACGTACGTAC";
        let mut g = Poa::new(seed, cfg());
        for _ in 0..4 {
            g.add_read(seed);
        }
        assert_eq!(g.consensus(), seed.to_vec());
    }

    #[test]
    fn majority_substitution_wins() {
        // 4 reads with C at pos 4, 1 read with A -> consensus keeps C
        let maj = b"ACGTCACGTA";
        let minr = b"ACGTAACGTA";
        let mut g = Poa::new(maj, cfg());
        g.add_read(maj);
        g.add_read(maj);
        g.add_read(minr);
        assert_eq!(g.consensus(), maj.to_vec());
    }

    #[test]
    fn read_paths_are_consistent_with_coverage() {
        // Provenance invariant: every node's `cov` equals the number of read
        // paths that visit it, and every path index is a valid node.
        let base = b"ACGTACGTAC";
        let ins = b"ACGTAACGTAC";
        let mut g = Poa::new(base, cfg());
        g.add_read(base);
        g.add_read(ins);
        g.add_read(base);
        assert_eq!(g.read_paths.len(), g.n_reads);
        let mut visit = vec![0u32; g.nodes.len()];
        for p in &g.read_paths {
            for &nd in p {
                assert!(nd < g.nodes.len());
                visit[nd] += 1;
            }
        }
        for (nd, &v) in visit.iter().enumerate() {
            assert_eq!(v, g.nodes[nd].cov, "node {nd} cov mismatch");
        }
    }

    #[test]
    fn insertion_by_minority_is_dropped() {
        let base = b"ACGTACGTAC";
        let ins = b"ACGTAACGTAC"; // extra A
        let mut g = Poa::new(base, cfg());
        g.add_read(base);
        g.add_read(base);
        g.add_read(ins);
        // heaviest path should follow the 3-read backbone, not the 1-read insert
        assert_eq!(g.consensus(), base.to_vec());
    }
}

#[cfg(test)]
mod tests_mismatch {
    use super::*;
    #[test]
    fn majority_substitution_over_seed_base() {
        // seed has A at pos 4; THREE reads have C there. Correct consensus = C.
        // If mismatches collapse onto the seed node, we'd wrongly keep A.
        let seed = b"ACGTAACGTA";
        let maj = b"ACGTCACGTA";
        let mut g = Poa::new(seed, PoaConfig::default());
        g.add_read(maj);
        g.add_read(maj);
        g.add_read(maj);
        let cons = g.consensus();
        assert_eq!(
            cons,
            maj.to_vec(),
            "expected majority C at pos4, got {}",
            String::from_utf8_lossy(&cons)
        );
    }
}

#[cfg(test)]
mod tests_semiglobal {
    use super::*;
    fn cfg() -> PoaConfig {
        PoaConfig::default()
    } // default = SemiGlobal

    #[test]
    fn partial_reads_align_to_their_subregion() {
        // Full locus; 3 full reads + 2 PARTIAL reads (left half, right half).
        // Semi-global should place partials in their subregion, not penalize the
        // uncovered flank, and the consensus should stay the full locus.
        let full = b"ACGTACGTTTGGCCAATTGGCCAAGTACAGTAC";
        let mut g = Poa::new(full, cfg());
        g.add_read(full);
        g.add_read(full);
        let left = &full[..16]; // covers only the 5' half
        let right = &full[16..]; // covers only the 3' half
        g.add_read(left);
        g.add_read(right);
        assert_eq!(
            g.consensus(),
            full.to_vec(),
            "partial reads must not truncate the consensus; got {}",
            String::from_utf8_lossy(&g.consensus())
        );
    }
}
