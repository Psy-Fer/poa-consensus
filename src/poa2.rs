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

/// Internal per-phase profiling (topo / align / integrate), gated behind the
/// `profile` feature so it is entirely absent — zero overhead — in normal builds.
/// The `profile`-on version accumulates elapsed nanoseconds per phase into atomics
/// the speed harness reads; the `profile`-off version is a no-op ZST guard.
/// Not part of the public API.
#[cfg(feature = "profile")]
#[doc(hidden)]
pub mod profile {
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::time::Instant;

    pub static TOPO_NS: AtomicU64 = AtomicU64::new(0);
    pub static ALIGN_NS: AtomicU64 = AtomicU64::new(0);
    pub static INTEGRATE_NS: AtomicU64 = AtomicU64::new(0);

    #[derive(Clone, Copy)]
    pub enum Phase {
        Topo,
        Align,
        Integrate,
    }

    pub struct Timer {
        start: Instant,
        phase: Phase,
    }
    pub fn guard(phase: Phase) -> Timer {
        Timer {
            start: Instant::now(),
            phase,
        }
    }
    impl Drop for Timer {
        fn drop(&mut self) {
            let ns = self.start.elapsed().as_nanos() as u64;
            let slot = match self.phase {
                Phase::Topo => &TOPO_NS,
                Phase::Align => &ALIGN_NS,
                Phase::Integrate => &INTEGRATE_NS,
            };
            slot.fetch_add(ns, Ordering::Relaxed);
        }
    }

    /// Reset all phase accumulators to zero.
    pub fn reset() {
        TOPO_NS.store(0, Ordering::SeqCst);
        ALIGN_NS.store(0, Ordering::SeqCst);
        INTEGRATE_NS.store(0, Ordering::SeqCst);
    }
    /// (topo_ns, align_ns, integrate_ns) accumulated since the last `reset`.
    pub fn take() -> (u64, u64, u64) {
        (
            TOPO_NS.load(Ordering::SeqCst),
            ALIGN_NS.load(Ordering::SeqCst),
            INTEGRATE_NS.load(Ordering::SeqCst),
        )
    }
}

#[cfg(not(feature = "profile"))]
mod profile {
    #[derive(Clone, Copy)]
    pub enum Phase {
        Topo,
        Align,
        Integrate,
    }
    /// No-op guard (ZST, no Drop) — compiles away entirely.
    pub struct NoOp;
    #[inline(always)]
    pub fn guard(_phase: Phase) -> NoOp {
        NoOp
    }
}

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
pub(crate) enum Op {
    Match(usize), // read base aligned to this node (match or mismatch)
    Ins(u8),      // read base with no graph node
    Del(usize),   // graph node skipped by the read
}

pub struct Poa {
    nodes: Vec<Node>,
    cfg: PoaConfig,
    pub(crate) n_reads: usize,
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
    pub(crate) fn topo_order(&self) -> (Vec<usize>, Vec<usize>) {
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
    pub(crate) fn align(&self, read: &[u8], topo: &[usize], rank_of: &[usize]) -> Vec<Op> {
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

        // Per-row DP scratch, reused across nodes (cleared + resized to each row's
        // `width` below). Hoisted out of the loop so the hot path does not allocate
        // 5 Vecs per node — the dominant allocation churn (see the `profile` phase
        // attribution: `align` is ~96% of build time).
        let mut mrow: Vec<i32> = Vec::new();
        let mut irow: Vec<i32> = Vec::new();
        let mut drow: Vec<i32> = Vec::new();
        let mut mbrow: Vec<(u8, u32)> = Vec::new();
        let mut dbrow: Vec<(u8, u32)> = Vec::new();

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

            mrow.clear();
            mrow.resize(width, NEG_INF);
            irow.clear();
            irow.resize(width, NEG_INF);
            drow.clear();
            drow.resize(width, NEG_INF);
            mbrow.clear();
            mbrow.resize(width, (0u8, VSTART));
            dbrow.clear();
            dbrow.resize(width, (0u8, VSTART));
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
                        // The forward recurrence is I[j] = max(M[j-1]+go+ge,
                        // I[j-1]+ge), so arriving from M pays the gap-open. The
                        // predecessor is M only when M[j] + go >= I[j]; comparing
                        // raw `o >= e` (go omitted) picks M too eagerly and cuts a
                        // multi-base insert run one base short, mis-emitting an
                        // inserted base as a Match. Guard the NEG_INF M cell so an
                        // unreachable M can never win the comparison.
                        state = if o != NEG_INF && o + go >= e { 0 } else { 1 };
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
        let (topo, rank_of) = {
            let _g = profile::guard(profile::Phase::Topo);
            self.topo_order()
        };
        let ops = {
            let _g = profile::guard(profile::Phase::Align);
            self.align(read, &topo, &rank_of)
        };
        let (path, medges) = {
            let _g = profile::guard(profile::Phase::Integrate);
            self.integrate(read, &ops)
        };
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
        // `visited` is cycle-safety insurance: the graph is a DAG by construction
        // (topo_order debug-asserts it), so `nxt` cannot cycle — but node fusion has
        // no explicit back-edge guard, so if a cycle ever slipped through in a
        // release build this bounds the walk instead of looping.
        let mut visited = vec![false; n];
        while cur != usize::MAX && !visited[cur] {
            visited[cur] = true;
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

    /// Fit score for a candidate consensus: build a throwaway graph seeded on
    /// `candidate`, align every read to it, and return the mean per-base total of
    /// Insert + Delete ops (content the candidate doesn't explain / doesn't
    /// confirm). Lower = better fit; 0.0 = every read matches exactly. Used to
    /// pick between heaviest-path and majority-frequency candidates per call.
    fn fit_score(candidate: &[u8], reads: &[&[u8]], cfg: &PoaConfig) -> f64 {
        if candidate.is_empty() {
            return f64::MAX;
        }
        let g = Poa::new(candidate, cfg.clone());
        let (topo, rank) = g.topo_order();
        let mut indel = 0usize;
        let mut bases = 0usize;
        for r in reads {
            for op in g.align(r, &topo, &rank) {
                if matches!(op, Op::Ins(_) | Op::Del(_)) {
                    indel += 1;
                }
            }
            bases += r.len();
        }
        indel as f64 / bases.max(1) as f64
    }

    /// Best of heaviest-path and majority-frequency, chosen by `fit_score`
    /// against the reads. The two consensus modes win different cases (heaviest
    /// on clean/short, majority on high-error length-variable repeats); scoring
    /// both against the actual reads and keeping the better fit gets the best of
    /// each without an oracle.
    pub fn consensus_best_fit(&self, reads: &[&[u8]]) -> Vec<u8> {
        let hp = self.consensus();
        let mf = self.consensus_majority();
        if hp == mf {
            return hp;
        }
        let fh = Self::fit_score(&hp, reads, &self.cfg);
        let fm = Self::fit_score(&mf, reads, &self.cfg);
        if fm < fh { mf } else { hp }
    }

    /// Rich `Consensus` output using the best-fit sequence (heaviest-path vs
    /// majority-frequency, chosen by `fit_score`). When heaviest-path wins (the
    /// common case) this is exactly `consensus_full()`. When majority-frequency
    /// wins, the emitted sequence is the MF one and the per-position fields
    /// (coverage, path_weights, gaps) are recomputed from the MF columns so the
    /// output stays self-consistent; graph-structural fields (graph_stats,
    /// bubble_sites) are unchanged since they describe the graph, not the path.
    pub fn consensus_full_best_fit(&self, reads: &[&[u8]]) -> crate::types::Consensus {
        let base = self.consensus_full();
        let (mf_seq, mf_cov) = self.consensus_majority_cov();
        if base.sequence == mf_seq {
            return base;
        }
        let fh = Self::fit_score(&base.sequence, reads, &self.cfg);
        let fm = Self::fit_score(&mf_seq, reads, &self.cfg);
        if fh <= fm {
            return base; // heaviest-path wins (ties keep the path-based output)
        }
        Self::mf_full(base, mf_seq, mf_cov)
    }

    /// Rich `Consensus` forced to the majority-frequency (MSA-column) sequence,
    /// regardless of fit. This is what `ConsensusMode::MajorityFrequency`
    /// selects — useful for high-depth amplicons where column majority is
    /// trusted over the heaviest path. Per-position fields come from the MF
    /// columns; graph-structural fields (graph_stats, bubble_sites) from the
    /// underlying graph.
    pub fn consensus_full_majority(&self) -> crate::types::Consensus {
        let base = self.consensus_full();
        let (mf_seq, mf_cov) = self.consensus_majority_cov();
        if base.sequence == mf_seq {
            return base;
        }
        Self::mf_full(base, mf_seq, mf_cov)
    }

    /// Rebuild a rich `Consensus` around a majority-frequency sequence `mf_seq`
    /// with per-column coverage `mf_cov`, taking graph-structural fields from
    /// `base` (the heaviest-path output over the same graph). Per-position
    /// fields are recomputed so the output stays self-consistent with the MF
    /// sequence: `path_weights` proxies incoming support as the min of adjacent
    /// column coverages (an MF consensus has no single node path), and interior
    /// coverage gaps are re-derived from `mf_cov`.
    fn mf_full(
        base: crate::types::Consensus,
        mf_seq: Vec<u8>,
        mf_cov: Vec<u32>,
    ) -> crate::types::Consensus {
        use crate::types::{CoverageGap, GapKind};
        let mut path_weights = Vec::with_capacity(mf_cov.len());
        for i in 0..mf_cov.len() {
            let w = if i == 0 {
                mf_cov[0]
            } else {
                mf_cov[i - 1].min(mf_cov[i])
            };
            path_weights.push(w as i32);
        }
        let mut gaps = Vec::new();
        let cl = mf_cov.len();
        let mut i = 0;
        while i < cl {
            if mf_cov[i] <= 1 {
                let s = i;
                while i < cl && mf_cov[i] <= 1 {
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
        crate::types::Consensus {
            sequence: mf_seq,
            coverage: mf_cov,
            path_weights,
            n_reads: base.n_reads,
            graph_stats: base.graph_stats,
            gaps,
            bubble_sites: base.bubble_sites,
            read_indices: vec![],
        }
    }

    /// Assign every node an MSA column index (union `aligned` siblings into one
    /// column, then topologically order the resulting column DAG). Returns
    /// (node_column, n_columns) or None if the column DAG is not acyclic (falls
    /// back to heaviest-path in that case).
    fn msa_columns(&self) -> Option<(Vec<usize>, usize)> {
        let n = self.nodes.len();
        // Union-find over aligned siblings.
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        for i in 0..n {
            for &a in &self.nodes[i].aligned {
                let (ri, ra) = (find(&mut parent, i), find(&mut parent, a));
                if ri != ra {
                    parent[ri] = ra;
                }
            }
        }
        // Group DAG: edges between column-groups (dedup, drop self-edges).
        let mut group_out: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut indeg = vec![0usize; n];
        let mut seen = std::collections::HashSet::new();
        for i in 0..n {
            let gi = find(&mut parent, i);
            for e in &self.nodes[i].out {
                let gj = find(&mut parent, e.to);
                if gi != gj && seen.insert((gi, gj)) {
                    group_out[gi].push(gj);
                    indeg[gj] += 1;
                }
            }
        }
        // Kahn on the group roots (a group's representative is its find() root).
        let roots: Vec<usize> = (0..n).filter(|&i| find(&mut parent, i) == i).collect();
        let mut stack: Vec<usize> = roots.iter().copied().filter(|&g| indeg[g] == 0).collect();
        stack.sort_unstable();
        let mut col_of_group = vec![usize::MAX; n];
        let mut ncol = 0usize;
        while let Some(g) = stack.pop() {
            col_of_group[g] = ncol;
            ncol += 1;
            let mut outs = group_out[g].clone();
            outs.sort_unstable();
            for gj in outs {
                indeg[gj] -= 1;
                if indeg[gj] == 0 {
                    stack.push(gj);
                }
            }
        }
        if ncol != roots.len() {
            return None; // cycle in the column DAG -> unusable
        }
        let mut node_col = vec![0usize; n];
        for (i, nc) in node_col.iter_mut().enumerate() {
            let g = find(&mut parent, i);
            *nc = col_of_group[g];
        }
        Some((node_col, ncol))
    }

    /// Most-frequent-base (column-majority MSA) consensus — SPOA/abPOA style.
    /// Each read votes its base at every column it occupies and a gap at every
    /// column it skips; a column emits its plurality base only if a real base
    /// beats the gap count. This counts deletions explicitly (unlike edge-weight
    /// heaviest-path), so a minority insertion becomes a gap-majority column and
    /// is correctly dropped -- the mechanism that fixes homopolymer over-call.
    pub fn consensus_majority(&self) -> Vec<u8> {
        self.consensus_majority_cov().0
    }

    /// Majority-frequency consensus plus per-emitted-column coverage (the count
    /// of reads present at that column). Falls back to heaviest-path (with its
    /// path coverage) if the column DAG is unusable.
    fn consensus_majority_cov(&self) -> (Vec<u8>, Vec<u32>) {
        let (node_col, ncol) = match self.msa_columns() {
            Some(x) => x,
            None => {
                let path = self.heaviest_path_nodes();
                let seq = path.iter().map(|&nd| self.nodes[nd].base).collect();
                let cov = path.iter().map(|&nd| self.nodes[nd].cov).collect();
                return (seq, cov);
            }
        };
        // Per-column base tallies (A,C,G,T) plus present-read count.
        let mut tally = vec![[0u32; 4]; ncol];
        let mut present = vec![0u32; ncol];
        // ACGT → bucket index; any other base (N / IUPAC) → None so it is NOT
        // folded into the T tally (it still counts as column occupancy below).
        let code = |b: u8| match b {
            b'A' => Some(0usize),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        };
        // `covering[c]` = reads whose occupied column *range* brackets column `c`
        // (via a difference array over each read's [min,max] column). A read that
        // simply doesn't reach `c` (a partial / non-spanning read) is NOT counted,
        // so under SemiGlobal it can't masquerade as a deletion and truncate the
        // consensus near the flanks. Only reads that span `c` but don't occupy it
        // are genuine deletions (`covering - present`).
        let mut cover_delta = vec![0i32; ncol + 1];
        for path in &self.read_paths {
            // A read occupies each column at most once; guard duplicates.
            let mut last_col = usize::MAX;
            let (mut lo, mut hi) = (usize::MAX, 0usize);
            for &nd in path {
                let c = node_col[nd];
                lo = lo.min(c);
                hi = hi.max(c);
                if c == last_col {
                    continue;
                }
                last_col = c;
                present[c] += 1; // column occupancy (any base, incl. N)
                if let Some(ci) = code(self.nodes[nd].base) {
                    tally[c][ci] += 1; // base plurality: ACGT only
                }
            }
            if lo != usize::MAX {
                cover_delta[lo] += 1;
                cover_delta[hi + 1] -= 1;
            }
        }
        let mut covering = vec![0u32; ncol];
        let mut acc = 0i32;
        for c in 0..ncol {
            acc += cover_delta[c];
            covering[c] = acc.max(0) as u32;
        }
        let bases = *b"ACGT";
        let mut out = Vec::with_capacity(ncol);
        let mut cov = Vec::with_capacity(ncol);
        // A column is emitted only if a majority of the whole read set REACHES it
        // (`covering >= floor`) AND its plurality base beats genuine deletions
        // among the reads that span it (`bc > gap`). The `covering` floor keeps
        // the majority-*length* semantics (a minority of longer reads can't emit
        // their private tail), while measuring `gap` against `covering` rather than
        // `n_reads` stops reads that simply END elsewhere from voting as deletions
        // and truncating a well-covered column.
        let floor = coverage_threshold(self.n_reads, self.cfg.min_coverage_fraction);
        for c in 0..ncol {
            let gap = covering[c] - present[c].min(covering[c]);
            let (bi, &bc) = tally[c]
                .iter()
                .enumerate()
                .max_by_key(|&(_, &v)| v)
                .unwrap();
            if bc > gap && covering[c] >= floor {
                out.push(bases[bi]);
                cov.push(bc); // reads supporting the emitted base at this column
            }
        }
        (out, cov)
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
    fn compute_stats(&self, er: &HashMap<(usize, usize), Vec<u32>>) -> crate::types::GraphStats {
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
                // Gini from the sorted-ascending Lorenz curve:
                //   G = (n+1)/n - 2*area/(n*T)  =  1 - 2*area/(n*T) + 1/n
                // The `+ 1/n` term was previously `- 1/n`, which biased every
                // graph by -2/n (a uniform graph scored -2/n instead of 0 and
                // could report negative inequality).
                1.0 - (2.0 * area as f64) / (nlen * tot as f64) + 1.0 / nlen
            }
        };
        // bubbles: fork nodes with 2+ out-edges each >= allele floor. Uses the
        // MATCHED edge view (`er`), not `Edge.weight`, so a length variant's
        // delete-bypass resume edge does not inflate bubble_count/max_bubble_depth
        // (see the matched-edge invariant and the bubble scan in `consensus_full`).
        let floor = ((self.n_reads as f64) * self.cfg.min_allele_freq).floor() as i32;
        let floor = floor.max(1);
        let mut bubble_count = 0usize;
        let mut max_bubble_depth = 0usize;
        let mut longest_span = 0usize;
        for (ni, nd) in self.nodes.iter().enumerate() {
            let mw = |to: usize| er.get(&(ni, to)).map_or(0, |v| v.len()) as i32;
            let sig: Vec<i32> = nd
                .out
                .iter()
                .map(|e| mw(e.to))
                .filter(|&w| w >= floor)
                .collect();
            if sig.len() >= 2 {
                bubble_count += 1;
                let mut s = sig.clone();
                s.sort_unstable_by(|a, b| b.cmp(a));
                max_bubble_depth = max_bubble_depth.max(s[1] as usize); // minority arm
                for e in &nd.out {
                    if mw(e.to) >= floor {
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
        // bubble sites on the consensus path. Arm significance/counts use the
        // MATCHED edge view (`read_matched_edges` via `edge_reads`), NOT
        // `Edge.weight`: a length variant's short-allele delete-bypass resume edge
        // carries `Edge.weight` but no matched reads, so counting `Edge.weight`
        // here would register it as a spurious competing arm. Single-allele
        // reporting; multi-allele phasing handles length variants via the
        // bypass-inclusive `read_paths` separately. See the `read_matched_edges`
        // invariant.
        let er = self.edge_reads();
        let floor = (((self.n_reads as f64) * self.cfg.min_allele_freq).floor() as i32).max(1);
        let mut bubble_sites = Vec::new();
        for (pos, &nd) in path.iter().enumerate() {
            let mw = |to: usize| er.get(&(nd, to)).map_or(0, |v| v.len()) as i32;
            let sig: Vec<&Edge> = self.nodes[nd]
                .out
                .iter()
                .filter(|e| mw(e.to) >= floor)
                .collect();
            if sig.len() >= 2 {
                let mut arm_read_counts = Vec::new();
                let mut arm_sequences = Vec::new();
                let mut is_structural = false;
                for e in &sig {
                    arm_read_counts.push(mw(e.to) as u32);
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
            graph_stats: self.compute_stats(&er),
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
    pub(crate) fn edge_reads(&self) -> HashMap<(usize, usize), Vec<u32>> {
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
    /// path takes, or MISSING if it doesn't traverse an out-edge of the fork or
    /// takes a below-floor minor branch.
    ///
    /// Uses `read_paths` (occupied-node adjacencies, INCLUDING delete-bypass
    /// "resume" edges), NOT the matched-only edges: a length variant's shorter
    /// allele diverges by *deleting* the tail and resuming — a bypass edge — so a
    /// matched-only view shows no fork there and misses the length split entirely.
    /// The bypass adjacency is exactly the shorter allele's arm.
    pub fn phasing_matrix(&self, floor: u32) -> crate::phasing::PhasingMatrix {
        use crate::phasing::{MISSING, PhasingMatrix, Site};
        // Per-read adjacency map (node -> successor) and edge weights, from the
        // occupied-node paths (bypass-inclusive).
        let mut read_next: Vec<std::collections::HashMap<usize, usize>> =
            vec![std::collections::HashMap::new(); self.n_reads];
        let mut mw: std::collections::HashMap<(usize, usize), u32> =
            std::collections::HashMap::new();
        for (r, path) in self.read_paths.iter().enumerate() {
            for w in path.windows(2) {
                read_next[r].insert(w[0], w[1]);
                *mw.entry((w[0], w[1])).or_default() += 1;
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
    pub(crate) fn matched_view(&self) -> (MatchedOut, MatchedInc) {
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
pub(crate) type MatchedOut = Vec<Vec<(usize, i32)>>;
/// Matched-only reverse adjacency: `inc[b] = [a, …]`.
pub(crate) type MatchedInc = Vec<Vec<usize>>;

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
    fn majority_mode_forces_msa_column_consensus() {
        // ConsensusMode::MajorityFrequency (via `consensus_full_majority`) emits
        // the MSA-column-majority sequence, not the best-fit/heaviest one. Even
        // when heaviest and majority agree here, the forced-majority rich output
        // must equal the raw column-majority sequence exactly (wiring guard for
        // the CLI `--consensus-mode majority` flag / config.consensus_mode).
        let seed = b"ACGTACGTACGTACGT";
        let ins = b"ACGTACGTAACGTACGT"; // one extra A (minority insertion)
        let mut g = Poa::new(seed, cfg());
        for _ in 0..3 {
            g.add_read(seed);
        }
        g.add_read(ins);
        assert_eq!(g.consensus_full_majority().sequence, g.consensus_majority());
        // The minority insertion is a gap-majority column, so it is dropped:
        // majority consensus equals the clean backbone.
        assert_eq!(g.consensus_majority(), seed.to_vec());
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

    #[test]
    #[ignore = "probe; run with --ignored --nocapture"]
    fn probe_homopolymer_dump() {
        // A homopolymer-bearing repeat (AAAAG = A4 G) where reads carry
        // per-unit A-run indel error, mimicking the AAAAG over-call. Dump the
        // graph so we can SEE whether insertion nodes proliferate into parallel
        // A-nodes vs fuse into a clean column chain.
        // Load the ACTUAL failing read set + truth from the compare grid.
        let dir = std::env::var("PROBE_DIR").unwrap_or_else(|_| "/tmp/poa_compare_grid".into());
        let name = std::env::var("PROBE_CASE").unwrap_or_else(|_| "AAAAG_u60_d20_e8_s1".into());
        let truth = std::fs::read(format!("{dir}/{name}.truth")).expect("truth file");
        let fa = std::fs::read_to_string(format!("{dir}/{name}.reads.fa")).expect("reads file");
        let mut reads: Vec<Vec<u8>> = Vec::new();
        let mut cur: Vec<u8> = Vec::new();
        for line in fa.lines() {
            if line.starts_with('>') {
                if !cur.is_empty() {
                    reads.push(std::mem::take(&mut cur));
                }
            } else {
                cur.extend_from_slice(line.trim().as_bytes());
            }
        }
        if !cur.is_empty() {
            reads.push(cur);
        }
        // Median-length seed, matching the CLI.
        let mut order: Vec<usize> = (0..reads.len()).collect();
        order.sort_by_key(|&i| reads[i].len());
        let seed_idx = order[order.len() / 2];
        let mut g = Poa::new(&reads[seed_idx], cfg());
        for (i, r) in reads.iter().enumerate() {
            if i != seed_idx {
                g.add_read(r);
            }
        }
        let cons = g.consensus();
        let ed = crate::poa2::tests::lev(&cons, &truth);
        let mf = g.consensus_majority();
        let edm = crate::poa2::tests::lev(&mf, &truth);
        let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let bf = g.consensus_best_fit(&refs);
        let edbf = crate::poa2::tests::lev(&bf, &truth);
        eprintln!(
            "case {name}: truth {} | heaviest {} edit {ed} | majority {} edit {edm} | bestfit edit {edbf}",
            truth.len(),
            cons.len(),
            mf.len()
        );
        // First divergence position (simple prefix scan).
        let div = cons
            .iter()
            .zip(truth.iter())
            .position(|(a, b)| a != b)
            .unwrap_or(0);
        eprintln!("first divergence at pos {div}");
        eprintln!(
            "  truth[{}..]: {}",
            div,
            String::from_utf8_lossy(&truth[div..(div + 40).min(truth.len())])
        );
        eprintln!(
            "  cons [{}..]: {}",
            div,
            String::from_utf8_lossy(&cons[div..(div + 40).min(cons.len())])
        );
        // Dump the heaviest-path nodes in a window around the divergence, with
        // their full out-edge fan (to see competing arms the picker chose among).
        let path = g.heaviest_path_nodes();
        let lo = div.saturating_sub(6);
        let hi = (div + 20).min(path.len());
        eprintln!("heaviest-path nodes [{lo}..{hi}] (base cov del | out-edges):");
        for (pi, &nd) in path.iter().enumerate().take(hi).skip(lo) {
            let n = &g.nodes[nd];
            let outs: Vec<String> = n
                .out
                .iter()
                .map(|e| format!("{}'{}'(w{})", e.to, g.nodes[e.to].base as char, e.weight))
                .collect();
            eprintln!(
                "  p{pi:<3} n{nd:<4} '{}' cov{:<2} del{:<2} -> {}",
                n.base as char,
                n.cov,
                n.del,
                outs.join(" ")
            );
        }
    }

    #[test]
    #[ignore = "grid sweep; run with --ignored --nocapture (one process, avoids per-exec sandbox tax)"]
    fn grid_sweep_consensus_modes() {
        let dir = std::env::var("PROBE_DIR").unwrap_or_else(|_| "/tmp/poa_compare_grid".into());
        let manifest = std::fs::read_to_string(format!("{dir}/manifest.txt")).expect("manifest");
        let (mut th, mut tm, mut tb) = (0usize, 0usize, 0usize);
        let mut n = 0usize;
        for name in manifest.lines().filter(|l| !l.trim().is_empty()) {
            let truth = match std::fs::read(format!("{dir}/{name}.truth")) {
                Ok(t) => t,
                Err(_) => continue,
            };
            let fa = std::fs::read_to_string(format!("{dir}/{name}.reads.fa")).unwrap();
            let mut reads: Vec<Vec<u8>> = Vec::new();
            let mut cur: Vec<u8> = Vec::new();
            for line in fa.lines() {
                if line.starts_with('>') {
                    if !cur.is_empty() {
                        reads.push(std::mem::take(&mut cur));
                    }
                } else {
                    cur.extend_from_slice(line.trim().as_bytes());
                }
            }
            if !cur.is_empty() {
                reads.push(cur);
            }
            if reads.len() < 3 {
                continue;
            }
            let mut order: Vec<usize> = (0..reads.len()).collect();
            order.sort_by_key(|&i| reads[i].len());
            let med = order[order.len() / 2];
            let mut g = Poa::new(&reads[med], cfg());
            for (i, r) in reads.iter().enumerate() {
                if i != med {
                    g.add_read(r);
                }
            }
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            th += lev(&g.consensus(), &truth);
            tm += lev(&g.consensus_majority(), &truth);
            tb += lev(&g.consensus_best_fit(&refs), &truth);
            n += 1;
        }
        eprintln!(
            "GRID SWEEP ({n} configs, full-Levenshtein vs truth):\n  heaviest={th}  majority={tm}  best_fit={tb}"
        );
    }

    // Op counts (ins, match) for aligning `read` against a graph seeded on `seed`.
    fn ins_match(seed: &[u8], read: &[u8]) -> (usize, usize) {
        let g = Poa::new(seed, cfg());
        let (topo, rank) = g.topo_order();
        let ops = g.align(read, &topo, &rank);
        let ins = ops.iter().filter(|o| matches!(o, Op::Ins(_))).count();
        let mat = ops.iter().filter(|o| matches!(o, Op::Match(_))).count();
        (ins, mat)
    }

    #[test]
    fn affine_insert_run_traceback_does_not_mislabel_inserts_as_matches() {
        // Regression for the affine insert-run traceback (forward recurrence pays
        // gap_open from M into I, so the predecessor is M only when M+go >= I;
        // the old `o >= e` test omitted go and picked M too eagerly, re-emitting
        // an inserted base as a Match).
        //
        // seed = "AAC", read = "ACCC": the correct alignment matches 2 nodes and
        // inserts 2 bases. The bug produced 3 matches / 1 insert — one inserted
        // 'C' laundered into a Match op, inflating a node's coverage. Verified
        // discriminating: this exact case flips (1,3)->(2,2) when the fix is
        // toggled (brute-force grid over all len-3/4 seeds x reads, 2026-07-29).
        let (ins, mat) = ins_match(b"AAC", b"ACCC");
        assert_eq!(
            (ins, mat),
            (2, 2),
            "inserted bases must stay Ins, not be laundered into Match ops"
        );
    }

    /// Small Levenshtein for the probe (no external dep).
    pub fn lev(a: &[u8], b: &[u8]) -> usize {
        let (n, m) = (a.len(), b.len());
        let mut prev: Vec<usize> = (0..=m).collect();
        let mut cur = vec![0usize; m + 1];
        for i in 1..=n {
            cur[0] = i;
            for j in 1..=m {
                let c = usize::from(a[i - 1] != b[j - 1]);
                cur[j] = (prev[j] + 1).min(cur[j - 1] + 1).min(prev[j - 1] + c);
            }
            std::mem::swap(&mut prev, &mut cur);
        }
        prev[m]
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
