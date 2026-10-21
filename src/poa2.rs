//! Partial-order alignment consensus engine (abPOA-style).
//!
//! Full partial-order affine-gap DP (one DP row per graph node in topological
//! order, every node relaxes ALL its in-edges; columns are read positions),
//! with an adaptive static-diagonal-union band and semi-global free ends by
//! default. Banded rows store only their window, 12 bytes per cell (M/I/D
//! `i32`); no backpointers are kept, the traceback replays the comparisons.
//!
//! Read integration uses SPOA/abPOA-style **node fusion** (match → reuse node;
//! mismatch → reuse an aligned sibling carrying that base, else create a new
//! substitution node and cross-link it into the alignment column; insert → new
//! node; delete → skip the node and count it in its `del`). Node fusion is in
//! the core deliberately: full DP does NOT obviate it (DP decides the
//! alignment; the graph update must still represent a mismatch as a
//! substitution node, or the read's base is dropped and the seed's minority
//! base is wrongly kept; see
//! `tests_mismatch::majority_substitution_over_seed_base`).
//!
//! Consensus: heaviest path, majority-frequency (MSA column plurality), or
//! best-fit between the two (`ConsensusMode`). Per-read `read_paths` /
//! `read_matched_edges` provenance feeds multi-allele phasing. The value the
//! crate adds lives ABOVE this engine (analysis, diagnostics, multi-allele);
//! the engine itself aims to match abPOA quality.
//!
//! Deferred: two-piece convex gaps and SIMD row scoring.

use crate::config::{AlignmentMode, PoaConfig};
use std::cell::Cell;
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

/// Highest support quality reported (Phred). Beyond this the binomial tail
/// says nothing useful: reads agreeing does not make them right.
pub(crate) const MAX_SUPPORT_QUALITY: u8 = 60;

/// Phred-scaled probability that a base's strongest alternative, not the
/// called base, is what the reads favour: with `agree` reads for the call and
/// `alt` for the alternative, each read a vote under a flat prior,
/// `P = P(Binomial(agree + alt + 1, 1/2) > agree)`. 10 agreeing reads and no
/// alternative give Q33; a 6-4 split gives Q6; a tie gives Q3.
fn support_phred(agree: u32, alt: u32) -> u8 {
    let n = (agree + alt + 1) as usize;
    let k_min = agree as usize + 1;
    let ln_2n = n as f64 * std::f64::consts::LN_2;
    // Sum C(n, k) / 2^n for k = n down to k_min, with ln C(n, k) updated by
    // C(n, k - 1) = C(n, k) * k / (n - k + 1).
    let (mut ln_c, mut p) = (0.0f64, 0.0f64);
    let mut k = n;
    loop {
        p += (ln_c - ln_2n).exp();
        if k == k_min {
            break;
        }
        ln_c += (k as f64).ln() - ((n - k + 1) as f64).ln();
        k -= 1;
    }
    let q = -10.0 * p.log10();
    q.clamp(0.0, MAX_SUPPORT_QUALITY as f64).round() as u8
}

/// Support quality for a consensus `seq` that is not a node path of the
/// graph it came from (a majority-frequency call, or a homopolymer splice):
/// build a support graph seeded with `seq`, add every read, and score the seed
/// path (the seed itself does not vote).
pub(crate) fn support_quality_by_realignment(
    seq: &[u8],
    reads: &[&[u8]],
    cfg: &PoaConfig,
) -> Vec<u8> {
    if seq.is_empty() {
        return Vec::new();
    }
    let cfg = PoaConfig {
        homopolymer_min_run: 0,
        homopolymer_gap: None,
        ..cfg.clone()
    };
    let mut g = Poa::new(seq, cfg);
    for r in reads {
        g.add_read(r);
    }
    let path: Vec<usize> = (0..seq.len()).collect();
    g.path_support_quality(&path, true)
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

/// Indel totals of all reads aligned to a candidate consensus (see
/// [`Poa::fit`]). Both candidates are scored against the same reads, so raw
/// totals compare directly.
#[derive(Clone, Copy)]
struct Fit {
    /// Sum over indel events of `|gap_open|` plus `|gap_extend|` per base.
    affine: u64,
    /// Indel bases (inserted plus deleted).
    indel_bases: u64,
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
    /// `align`'s banded M/I/D matrices, kept between calls so their capacity
    /// is reused: the band is wider than its estimate wherever the adaptive
    /// and static diagonals diverge, and regrowing per read was ~9% of runtime.
    align_bufs: Cell<AlignBufs>,
    /// Gap costs `(open, extend)` for an indel that copies its neighbouring
    /// unit in a locally periodic stretch of the read (periods 2-6), or `None`.
    /// Set only on candidate builds (`with_tandem_gap`).
    tandem_gap: Option<(i32, i32)>,
}

#[derive(Default)]
struct AlignBufs {
    m: Vec<i32>,
    i: Vec<i32>,
    d: Vec<i32>,
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
            align_bufs: Cell::default(),
            tandem_gap: None,
        }
    }

    /// Price indels that copy a neighbouring repeat unit at `(open, extend)`
    /// (see the tandem context in `align`).
    pub(crate) fn with_tandem_gap(mut self, costs: (i32, i32)) -> Self {
        self.tandem_gap = Some(costs);
        self
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

    /// Full partial-order affine-gap alignment of `read` to the graph, global
    /// or semi-global per `alignment_mode`, banded per the band settings.
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
        // Homopolymer-aware gap costs (`PoaConfig::homopolymer_gap`): an indel
        // inside a homopolymer pays (hgo, hge) instead of (go, ge). Both sides
        // must really be a run, or a periodic repeat's unit indel (CAG|CAG, the
        // deleted C next to the read's next C) would count and get cheap.
        // Insertion of read[j-1]: that base is part of a run in the read.
        // Deletion of node t at read column j: t has a same-base graph
        // neighbour, and the read has that base adjacent (read[j-1] or read[j]).
        let hp = self.cfg.homopolymer_gap;
        let ins_run: Vec<bool> = (0..=l)
            .map(|j| {
                hp.is_some()
                    && j >= 1
                    && ((j >= 2 && read[j - 1] == read[j - 2]) || (j < l && read[j - 1] == read[j]))
            })
            .collect();
        let node_run: Vec<bool> = if hp.is_some() {
            self.nodes
                .iter()
                .map(|nd| {
                    nd.inc.iter().any(|&p| self.nodes[p].base == nd.base)
                        || nd.out.iter().any(|e| self.nodes[e.to].base == nd.base)
                })
                .collect()
        } else {
            Vec::new()
        };
        // Per-column insert costs, and the read padded by one sentinel on each
        // side (`rpad[j] = read[j-1]`, outside 1..=l a value no u8 base equals)
        // so the fill's delete-cost test needs no bounds checks.
        let tg = self.tandem_gap;
        // Periodic stretches of the read for periods 2..=6: wherever a unit is
        // immediately followed by a copy of itself, both copies are flagged, and
        // each read column inside such a stretch records which bases its unit
        // contains (bit per A/C/G/T). Purely local: no tract or unit is
        // identified, only "an indel here would copy its neighbour".
        let base_bit = |b: u8| match b {
            b'A' | b'a' => 1u8,
            b'C' | b'c' => 2,
            b'G' | b'g' => 4,
            b'T' | b't' => 8,
            _ => 0,
        };
        let mut tpos = vec![false; l];
        let mut tmask = vec![0u8; l + 1]; // per column j (between read[j-1] and read[j])
        if tg.is_some() {
            for per in 2..=6usize {
                if l < 2 * per {
                    break;
                }
                for st in 0..=l - 2 * per {
                    if read[st..st + per] == read[st + per..st + 2 * per] {
                        let m = read[st..st + per].iter().fold(0u8, |m, &b| m | base_bit(b));
                        tpos[st..st + 2 * per].fill(true);
                        // interior columns of the stretch
                        for x in &mut tmask[st + 1..st + 2 * per] {
                            *x |= m;
                        }
                    }
                }
            }
        }
        let icost: Vec<(i32, i32)> = ins_run
            .iter()
            .enumerate()
            .map(|(j, &run)| match (hp, tg) {
                (Some(c), _) if run => c,
                (_, Some(c)) if j >= 1 && tpos[j - 1] => c,
                _ => (go, ge),
            })
            .collect();
        let ins_cost = |j: usize| icost[j];
        const NO_BASE: u16 = 0x100;
        let rpad: Vec<u16> = std::iter::once(NO_BASE)
            .chain(read.iter().map(|&b| b as u16))
            .chain(std::iter::once(NO_BASE))
            .collect();
        let del_cost = |t: usize, j: usize, b: u8| match (hp, tg) {
            (Some(c), _) if node_run[t] && (rpad[j] == b as u16 || rpad[j + 1] == b as u16) => c,
            (_, Some(c)) if tmask[j] & base_bit(b) != 0 => c,
            _ => (go, ge),
        };

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
        // Half-width: `band_width`, widened to the abPOA length term
        // `b + ceil(f * L)` when `adaptive_band` is on. With it off the band is
        // the fixed `band_width`.
        let half = {
            let adaptive = if self.cfg.adaptive_band {
                self.cfg.adaptive_band_b + (self.cfg.adaptive_band_f * l as f32).ceil() as usize
            } else {
                0
            };
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
        let AlignBufs {
            m: mut mflat,
            i: mut iflat,
            d: mut dflat,
        } = self.align_bufs.take();
        for v in [&mut mflat, &mut iflat, &mut dflat] {
            v.clear();
            v.reserve(cap);
        }
        let mut d0 = vec![NEG_INF; n];
        let mut d0bk = vec![(0u8, VSTART); n];
        let mut best_j = vec![1usize; n]; // per-rank winning column (adaptive center)

        // Per-row accumulators over the predecessors, reused across nodes:
        // `macc[k]` = max of a predecessor's M/I/D at column j-1 (feeds M[j]),
        // `oacc[k]` = max of its M/I at j (D gap-open), `eacc[k]` = max of its
        // D at j (D extend). Every stored cell is exactly NEG_INF or a real
        // score, so a raw max is NEG_INF iff no predecessor reaches the cell,
        // and adding the cell's (shared) cost afterwards gives the same value
        // as the per-candidate comparison. Predecessor-outer / column-inner
        // keeps these loops branch-free over contiguous slices.
        let mut macc: Vec<i32> = Vec::new();
        let mut oacc: Vec<i32> = Vec::new();
        let mut eacc: Vec<i32> = Vec::new();

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
                let (dgo, dge) = del_cost(t, 0, base);
                if is_source {
                    // leading deletion from the virtual start
                    if dgo + dge > best {
                        best = dgo + dge;
                        bk = (0, VSTART);
                    }
                }
                for &p in &self.nodes[t].inc {
                    let pr = rank_of[p] as u32;
                    let de = d0[rank_of[p]];
                    if de != NEG_INF && de + dge > best {
                        best = de + dge;
                        bk = (2, pr);
                    }
                }
                d0[r] = best;
                d0bk[r] = bk;
            }

            macc.clear();
            macc.resize(width, NEG_INF);
            oacc.clear();
            oacc.resize(width, NEG_INF);
            eacc.clear();
            eacc.resize(width, NEG_INF);
            for &p in &self.nodes[t].inc {
                let pr = rank_of[p];
                let (plo, phi, po) = (lo_[pr], hi_[pr], off[pr]);
                let pw = phi + 1 - plo;
                let (pm, pi, pd) = (
                    &mflat[po..po + pw],
                    &iflat[po..po + pw],
                    &dflat[po..po + pw],
                );
                // M[j] reads the predecessor at j-1; j-1 = 0 is the D0 column.
                if lo == 1 {
                    macc[0] = macc[0].max(d0[pr]);
                }
                let (a, b) = (lo.max(plo + 1), hi.min(phi + 1));
                if a <= b {
                    let (k0, p0, len) = (a - lo, a - 1 - plo, b + 1 - a);
                    let acc = &mut macc[k0..k0 + len];
                    let src = pm[p0..p0 + len]
                        .iter()
                        .zip(&pi[p0..p0 + len])
                        .zip(&pd[p0..p0 + len]);
                    for (x, ((&m, &i), &d)) in acc.iter_mut().zip(src) {
                        *x = (*x).max(m.max(i).max(d));
                    }
                }
                // D[j] reads the predecessor at the same j.
                let (a, b) = (lo.max(plo), hi.min(phi));
                if a <= b {
                    let (k0, p0, len) = (a - lo, a - plo, b + 1 - a);
                    let open = oacc[k0..k0 + len]
                        .iter_mut()
                        .zip(pm[p0..p0 + len].iter().zip(&pi[p0..p0 + len]));
                    for (x, (&m, &i)) in open {
                        *x = (*x).max(m.max(i));
                    }
                    for (x, &d) in eacc[k0..k0 + len].iter_mut().zip(&pd[p0..p0 + len]) {
                        *x = (*x).max(d);
                    }
                }
            }

            // ---- M[t][j] and D[t][j] ---- No dependency along the row, so this
            // pass is branch-free (selects only) and vectorises.
            let o = mflat.len();
            mflat.resize(o + width, NEG_INF);
            iflat.resize(o + width, NEG_INF);
            dflat.resize(o + width, NEG_INF);
            let (mrow, irow, drow) = (&mut mflat[o..], &mut iflat[o..], &mut dflat[o..]);
            // Delete costs: homopolymer-priced where t is in a graph run and the
            // read has t's base adjacent to column j (see `del_cost`).
            let hp_row = hp.filter(|_| node_run[t]);
            let (hgo, hge) = hp_row.unwrap_or((go, ge));
            let bb = base as u16;
            let tmask_row = &tmask[lo..=hi];
            let tbit = if tg.is_some() { base_bit(base) } else { 0 };
            let (tgo, tge) = tg.unwrap_or((go, ge));
            let cols = rpad[lo..=hi]
                .iter()
                .zip(&rpad[lo + 1..=hi + 1])
                .zip(tmask_row);
            let accs = macc.iter().zip(oacc.iter().zip(&eacc));
            let outs = mrow.iter_mut().zip(drow.iter_mut());
            for (((m, d), (&ma, (&oa, &ea))), ((&rb, &rn), &tm)) in outs.zip(accs).zip(cols) {
                let sc = if rb == bb { mm } else { mis };
                *m = if ma == NEG_INF { NEG_INF } else { ma + sc };
                let in_run = hp_row.is_some() && (rb == bb || rn == bb);
                let in_tandem = tm & tbit != 0;
                let (dgo, dge) = if in_run {
                    (hgo, hge)
                } else if in_tandem {
                    (tgo, tge)
                } else {
                    (go, ge)
                };
                let d_open = if oa == NEG_INF {
                    NEG_INF
                } else {
                    oa + dgo + dge
                };
                let d_ext = if ea == NEG_INF { NEG_INF } else { ea + dge };
                *d = d_open.max(d_ext);
            }
            if (is_source || semi) && lo == 1 {
                // free start: read[0] begins the alignment here
                let sc = if base == read[0] { mm } else { mis };
                mrow[0] = mrow[0].max(sc);
            }

            // ---- I[t][j] ---- (read base inserted; same row, j-1) and the
            // row's best column. Sequential along the row.
            let mut row_best = NEG_INF;
            let mut row_best_j = best_j[r];
            let (mut m_prev, mut i_prev) = (NEG_INF, NEG_INF);
            for (k, (&(igo, ige), (&m, (ins_out, &d)))) in icost[lo..=hi]
                .iter()
                .zip(mrow.iter().zip(irow.iter_mut().zip(drow.iter())))
                .enumerate()
            {
                let i_open = if m_prev == NEG_INF {
                    NEG_INF
                } else {
                    m_prev + igo + ige
                };
                let i_ext = if i_prev == NEG_INF {
                    NEG_INF
                } else {
                    i_prev + ige
                };
                let ins = i_open.max(i_ext);
                *ins_out = ins;
                let cell = m.max(ins).max(d);
                if cell != NEG_INF && cell > row_best {
                    row_best = cell;
                    row_best_j = lo + k;
                }
                (m_prev, i_prev) = (m, ins);
            }
            best_j[r] = row_best_j;
        }

        // Backpointers are not stored during the fill: traceback re-derives the
        // winning predecessor of an M / D cell by replaying the fill's
        // comparisons in the same order (free start, then each predecessor's
        // M, I, D, strict `>`), so ties resolve exactly as a per-cell record
        // would. This is O(in-degree) per traceback step instead of 16 bytes
        // written per DP cell.
        let m_back = |r: usize, j: usize| -> (u8, u32) {
            let t = topo[r];
            let sc = if self.nodes[t].base == read[j - 1] {
                mm
            } else {
                mis
            };
            let mut best = NEG_INF;
            let mut bk = (0u8, VSTART);
            if (self.nodes[t].inc.is_empty() || semi) && j == 1 {
                best = sc;
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
            bk
        };
        let d_back = |r: usize, j: usize| -> (u8, u32) {
            let t = topo[r];
            let (dgo, dge) = del_cost(t, j, self.nodes[t].base);
            let mut best = NEG_INF;
            let mut bk = (0u8, VSTART);
            for &p in &self.nodes[t].inc {
                let pr = rank_of[p];
                let mo = bget(&mflat, &off, &lo_, &hi_, pr, j);
                let io = bget(&iflat, &off, &lo_, &hi_, pr, j);
                let de = bgetd(&dflat, &d0, &off, &lo_, &hi_, pr, j);
                for (st, v, c) in [(0u8, mo, dgo + dge), (1u8, io, dgo + dge), (2u8, de, dge)] {
                    if v != NEG_INF && v + c > best {
                        best = v + c;
                        bk = (st, pr as u32);
                    }
                }
            }
            bk
        };

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
                    let (pst, pr) = m_back(r, j);
                    j -= 1;
                    state = pst;
                    rr = pr;
                }
                1 => {
                    ops.push(Op::Ins(read[j - 1]));
                    // Open cost of THIS insert cell (must match the fill).
                    let igo = ins_cost(j).0;
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
                        state = if o != NEG_INF && o + igo >= e { 0 } else { 1 };
                    }
                }
                _ => {
                    ops.push(Op::Del(t));
                    let (pst, pr) = if j == 0 { d0bk[r] } else { d_back(r, j) };
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
        self.align_bufs.set(AlignBufs {
            m: mflat,
            i: iflat,
            d: dflat,
        });
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

    /// Fit of a candidate consensus to the reads: build a throwaway graph
    /// seeded on `candidate`, align every read to it, and total the indels.
    /// Mismatches are not counted: on noisy reads they are mostly read error,
    /// and counting them let a frame-broken majority candidate win on ONT
    /// CAG50. Used to pick between heaviest-path and majority-frequency
    /// candidates (see [`Poa::mf_beats_hp`]).
    fn fit(candidate: &[u8], reads: &[&[u8]], cfg: &PoaConfig) -> Fit {
        if candidate.is_empty() {
            return Fit {
                affine: u64::MAX,
                indel_bases: u64::MAX,
            };
        }
        let g = Poa::new(candidate, cfg.clone());
        let (topo, rank) = g.topo_order();
        let (open, extend) = (
            cfg.gap_open.unsigned_abs() as u64,
            cfg.gap_extend.unsigned_abs() as u64,
        );
        let mut fit = Fit {
            affine: 0,
            indel_bases: 0,
        };
        for r in reads {
            let mut prev = 0u8; // 0 = match, 1 = insert, 2 = delete
            for op in g.align(r, &topo, &rank) {
                let kind = match op {
                    Op::Match(_) => 0,
                    Op::Ins(_) => 1,
                    Op::Del(_) => 2,
                };
                if kind != 0 {
                    fit.indel_bases += 1;
                    fit.affine += extend + if kind != prev { open } else { 0 };
                }
                prev = kind;
            }
        }
        fit
    }

    /// Whether the majority-frequency candidate should replace the heaviest
    /// path: only when it fits the reads strictly better on BOTH measures,
    /// fewer indel bases and a lower affine indel penalty. The heaviest path is
    /// the graph-supported default; MF overrides it only on unambiguous
    /// evidence. Each measure alone fails on (CA)n ladders:
    /// - Indel bases alone reward a candidate halfway between two repeat
    ///   lengths: with reads at 19 and 20 units, a frame-broken `(CA)19 C` is
    ///   one base from every read, so it beat both real lengths.
    /// - The affine penalty alone breaks exact base-count ties toward the
    ///   shorter candidate: a read longer than the candidate can absorb several
    ///   differences into one insertion run (one open), so a one-unit-short MF
    ///   won on the open charge alone (and a 6-C homopolymer became 5 C on
    ///   HG002 HiFi DM1).
    fn mf_beats_hp(hp: &[u8], mf: &[u8], reads: &[&[u8]], cfg: &PoaConfig) -> bool {
        let fh = Self::fit(hp, reads, cfg);
        let fm = Self::fit(mf, reads, cfg);
        fm.affine < fh.affine && fm.indel_bases < fh.indel_bases
    }

    /// Best of heaviest-path and majority-frequency, chosen by
    /// [`Poa::mf_beats_hp`]. The two consensus modes win different cases
    /// (heaviest on clean/short, majority on high-error length-variable
    /// repeats); scoring both against the actual reads gets the better of each
    /// without an oracle.
    pub fn consensus_best_fit(&self, reads: &[&[u8]]) -> Vec<u8> {
        let hp = self.consensus();
        let mf = self.consensus_majority();
        if hp != mf && Self::mf_beats_hp(&hp, &mf, reads, &self.cfg) {
            mf
        } else {
            hp
        }
    }

    /// Rich `Consensus` output using the best-fit sequence (heaviest-path vs
    /// majority-frequency, chosen by [`Poa::mf_beats_hp`]). When heaviest-path
    /// wins (the common case) this is exactly `consensus_full()`. When
    /// majority-frequency wins, the emitted sequence is the MF one and the
    /// per-position fields (coverage, path_weights, gaps) are recomputed from
    /// the MF columns so the output stays self-consistent; graph-structural
    /// fields (graph_stats, bubble_sites) are unchanged since they describe the
    /// graph, not the path.
    pub fn consensus_full_best_fit(&self, reads: &[&[u8]]) -> crate::types::Consensus {
        let base = self.consensus_full();
        let (mf_seq, mf_cov) = self.consensus_majority_cov();
        if base.sequence == mf_seq || !Self::mf_beats_hp(&base.sequence, &mf_seq, reads, &self.cfg)
        {
            return base;
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
            // Not a node path in this graph: the caller scores it against the
            // reads (`support_quality_by_realignment`).
            support_quality: Vec::new(),
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
        let (seq, cov, _) = self.consensus_majority_cols();
        (seq, cov)
    }

    /// `consensus_majority_cov` plus the MSA column of each emitted base
    /// (`None` when it fell back to the heaviest path).
    fn consensus_majority_cols(&self) -> (Vec<u8>, Vec<u32>, Option<Vec<usize>>) {
        let (node_col, ncol) = match self.msa_columns() {
            Some(x) => x,
            None => {
                let path = self.heaviest_path_nodes();
                let seq = path.iter().map(|&nd| self.nodes[nd].base).collect();
                let cov = path.iter().map(|&nd| self.nodes[nd].cov).collect();
                return (seq, cov, None);
            }
        };
        // Per-column symbol tallies (A,C,G,T, then the homopolymer-compressed
        // poly-A..poly-T symbols, which set the high bit; see `crate::hpc`)
        // plus present-read count.
        let mut tally = vec![[0u32; 8]; ncol];
        let mut present = vec![0u32; ncol];
        // ACGT (± poly bit) → bucket index; any other base (N / IUPAC) → None
        // so it is NOT folded into a real tally (it still counts as column
        // occupancy below).
        let code = |b: u8| {
            let i = match b & 0x7f {
                b'A' => 0usize,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
            Some(if b & 0x80 != 0 { i + 4 } else { i })
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
                    tally[c][ci] += 1; // base plurality: ACGT (+ poly) only
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
        let bases = [
            b'A',
            b'C',
            b'G',
            b'T',
            b'A' | 0x80,
            b'C' | 0x80,
            b'G' | 0x80,
            b'T' | 0x80,
        ];
        let mut out = Vec::with_capacity(ncol);
        let mut cov = Vec::with_capacity(ncol);
        let mut cols = Vec::with_capacity(ncol);
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
                cols.push(c);
            }
        }
        (out, cov, Some(cols))
    }

    /// For a consensus `cons` of this graph (its heaviest path or its
    /// majority-frequency sequence): the MSA column of every node, and the
    /// column of every consensus position. `None` if `cons` is neither or the
    /// column DAG is unusable.
    pub(crate) fn consensus_columns(&self, cons: &[u8]) -> Option<(Vec<usize>, Vec<usize>)> {
        let (node_col, _) = self.msa_columns()?;
        let path = self.heaviest_path_nodes();
        if path.len() == cons.len()
            && path
                .iter()
                .zip(cons)
                .all(|(&n, &b)| self.nodes[n].base == b)
        {
            let pos_col = path.iter().map(|&n| node_col[n]).collect();
            return Some((node_col, pos_col));
        }
        let (mf, _, cols) = self.consensus_majority_cols();
        if mf == cons {
            return Some((node_col, cols?));
        }
        None
    }

    /// Occupied nodes per read, in insertion order (seed first).
    pub(crate) fn read_paths(&self) -> &[Vec<usize>] {
        &self.read_paths
    }

    pub(crate) fn node_base(&self, n: usize) -> u8 {
        self.nodes[n].base
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
    /// Minimum matched-read count for a bubble arm to be reported (bubble
    /// sites, bubble statistics): `ceil(n_reads * min_allele_freq)`, at least
    /// 1, so a reported arm carries at least `min_allele_freq` of the reads.
    /// Multi-allele discovery uses the same rule; acting on a split
    /// additionally needs 2 reads (`multi::phasing_floor`).
    fn allele_arm_floor(&self) -> i32 {
        (((self.n_reads as f64) * self.cfg.min_allele_freq).ceil() as i32).max(1)
    }

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
        let floor = self.allele_arm_floor();
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
    /// Per-base support quality of a consensus that is the node path `path`
    /// in this graph (see [`Consensus::support_quality`]). Each read votes from
    /// its recorded path: between two consecutive consensus nodes it visits
    /// (`a`, then `b`), stepping straight from `a` to `a + 1` agrees with base
    /// `a + 1`; any other route is an alternative, against the skipped bases
    /// `a + 1 .. b` or, for an insertion (`b == a + 1` via other nodes), against
    /// `b`. Reads taking the same route (same anchors and intermediate nodes)
    /// form one alternative, and each base is scored against its strongest
    /// alternative only: sequencing errors scatter over many one-read routes,
    /// while a competing call concentrates on one. Bases outside a read's first
    /// and last consensus node get no vote from it. `skip_seed` drops read 0
    /// (a support graph seeded with the consensus itself).
    ///
    /// [`Consensus::support_quality`]: crate::types::Consensus::support_quality
    pub(crate) fn path_support_quality(&self, path: &[usize], skip_seed: bool) -> Vec<u8> {
        let n = path.len();
        let mut idx = vec![usize::MAX; self.nodes.len()];
        for (i, &nd) in path.iter().enumerate() {
            idx[nd] = i;
        }
        let mut agree = vec![0u32; n];
        let mut routes: HashMap<(usize, usize, Vec<usize>), u32> = HashMap::new();
        for (r, rp) in self.read_paths.iter().enumerate() {
            if skip_seed && r == 0 {
                continue;
            }
            // (consensus index, position in the read's path) of the last anchor
            let mut prev: Option<(usize, usize)> = None;
            for (k, &nd) in rp.iter().enumerate() {
                let ci = idx[nd];
                if ci == usize::MAX {
                    continue;
                }
                match prev {
                    None => agree[ci] += 1,
                    Some((a, _)) if ci <= a => continue, // not forward along the path
                    Some((a, ka)) => {
                        if ci == a + 1 && k == ka + 1 {
                            agree[ci] += 1;
                        } else {
                            *routes.entry((a, ci, rp[ka + 1..k].to_vec())).or_insert(0) += 1;
                            if ci > a + 1 {
                                // skipped or substituted bases are the alternative;
                                // the read still carries base `ci` itself
                                agree[ci] += 1;
                            }
                        }
                    }
                }
                prev = Some((ci, k));
            }
        }
        let mut alt = vec![0u32; n];
        for ((a, b, _), &count) in &routes {
            let span = if *b > a + 1 { a + 1..*b } else { *b..b + 1 };
            for x in &mut alt[span] {
                *x = (*x).max(count);
            }
        }
        agree
            .iter()
            .zip(&alt)
            .map(|(&a, &m)| support_phred(a, m))
            .collect()
    }

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
        let floor = self.allele_arm_floor();
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
            support_quality: self.path_support_quality(&path, false),
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
mod tests_best_fit_frame {
    use super::*;

    /// ONT (CA)20 reads from bladerunner / str-replicant (synthetic): four
    /// exact, three one CA short, the median-length seed (r4) missing one A.
    /// Majority-frequency splits the column vote where the short reads place
    /// their deletion and emits a frame-broken `(CA)19 C`; an indel-*bases*
    /// fit score preferred it (one base from every read). The heaviest path is
    /// right, and BestFit must keep it in every read order.
    const READS: [&[u8]; 10] = [
        b"GGGAAACTATTGTGCAATTTCACACACACACACACACACACACACACACACACACACACAGGATCCCATGACTCATCACG",
        b"GGGAAACTATTGTGCAATTTCACACACACACACACACACACACACACACACACACACACAGGATCCCATGACTCATCACG",
        b"GGGAAACTATTGTGCAATTTCACACACACACACACACACACACACACACACACACACAGGATCCCATGACTCATCACG",
        b"GGGAAACTATTGTGCAATTTCACACACACACACACACACACACCACACACACACACAGGATCCCATGACTCATCACG",
        b"GGGAAACTATTGTGCAATTTCACACACACACACACACACCACACACACACACACACACAGGATCCCATGACTCATCACG",
        b"GGGAAACTATTGTGCAATTCACACACACACACACACAACACACACACACACACACACAGGATCCATGACTCATCACG",
        b"GGGAAACTATTGTGCAATTTCACACACACACACACACACACACACACACACACACACAGGATCCCATGACTCATCACG",
        b"GGGAAACTATTGTGCAATTTCACACACACACACACACACACACACACACACACACACACAGGATCCCATGACTCATCACG",
        b"GGGAAACTATTGTGCAATTTCACACACACACACACACACACACACACACACACACACAGGATCCCATGACTCATCAGG",
        b"GGGAAACTATTGTGCAATTTCACACACACACACACACACACACACACACACACACACACAGGATCCCATGACTCATCACG",
    ];
    const TRUTH: &[u8] =
        b"GGGAAACTATTGTGCAATTTCACACACACACACACACACACACACACACACACACACACAGGATCCCATGACTCATCACG";

    #[test]
    fn best_fit_rejects_frame_broken_majority_in_every_order() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        for trial in 0..50 {
            let mut order: Vec<usize> = (0..READS.len()).collect();
            for i in (1..order.len()).rev() {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                order.swap(i, (state % (i as u64 + 1)) as usize);
            }
            let reads: Vec<&[u8]> = order.iter().map(|&i| READS[i]).collect();
            let c = crate::consensus(&reads, 0, &PoaConfig::default()).unwrap();
            assert_eq!(
                c.sequence,
                TRUTH,
                "trial {trial} order {order:?}: got {}",
                String::from_utf8_lossy(&c.sequence)
            );
        }
    }

    /// ONT (CA)22 reads from bladerunner / str-replicant (synthetic), split
    /// 4 x 21 units against 3 x 22 plus noise. Both candidates have the same
    /// indel *bases*; the one-unit-short majority candidate needs one fewer
    /// indel *event* (a longer read absorbs its differences into one insert
    /// run), so an affine-only fit picked it. MF must not win on a tie in
    /// indel bases.
    #[test]
    fn best_fit_keeps_heaviest_path_on_indel_base_tie() {
        let reads: [&[u8]; 10] = [
            b"CGGGTCGGTCTTTGCGTTTTCACACACACACACACACACACACACACACACACACACACACAGTTGGGAGTAGTGTTCTAGT",
            b"CGGGTCGGTCTTTGCGTTTTCACACACACACACACACACACACACACACACACACACACACACAGTTGGGAGTAGTGTTCTAG",
            b"CGGGTCGGTCTTTGCGTTTTCACACGCACACACACACACACACTCACACGCACACACACACAGTTGGGAGTAGTGTTCTATT",
            b"CGGGTCGGTCTTTGCGTTTTCACACACACACACACACACACACACACACACACACACACACAGTTGGGAGTAGTGTTCTAGT",
            b"CGGGTCGGTCTTTGCGTTTTCACAGACACACACACACACACACACAAGCACACACACACACACACAGTTGGGAGTAGTGTTCTAGT",
            b"CGGGTCGGTCTTTGCGTTTTCACACACACACACACACACACACACACACACACACACACACACAGTTGGGAGTAGTGTTCTAGT",
            b"CGGGTCGGTCTTTGCGTTTTCACACACACACACACACACACACACACACACACACACACACAGTTGGGAGTAGTGTTCTAGT",
            b"CGGGTCGGTCTTTGCGTTTTCACACACACACACACACACACACACACACACACACACACACAGTTGGGAGTAGTGTTCTAGT",
            b"CGGGTCGGTCTTTGCGTTTTCACACACACACACACACACACACACACACACACACACACACACACAGTTGGGAGTAGTGTTCTAGT",
            b"CGGGTCGGTCTTTGCGTTTTCACACACACACACACACACACACACACACACACACACACACACAGTTGGGAGTAGTGTTCTAGT",
        ];
        let truth: &[u8] =
            b"CGGGTCGGTCTTTGCGTTTTCACACACACACACACACACACACACACACACACACACACACACAGTTGGGAGTAGTGTTCTAGT";
        let c = crate::consensus(&reads, 0, &PoaConfig::default()).unwrap();
        assert_eq!(
            c.sequence,
            truth,
            "got {}",
            String::from_utf8_lossy(&c.sequence)
        );
    }

    /// Synthetic analogue of HG002 HiFi DM1 (homozygous `C6 (CAG)11`): two
    /// reads one C short, one exact, and one exact plus an extra CAG. Against
    /// a 5-C candidate the long read takes its extra C and CAG as one adjacent
    /// 4-base insertion, so the 5-C candidate saves an indel open at the same
    /// indel-base count. The heaviest path (6 C) must be kept.
    #[test]
    fn best_fit_keeps_longer_homopolymer_on_indel_base_tie() {
        let reads: [&[u8]; 4] = [
            b"GCTAAAGACAATTACATAATCCCCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTACACGTCAGCACGAAACT",
            b"GCTAAAGACAATTACATAATCCCCCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTACACGTCAGCACGAAACT",
            b"GCTAAAGACAATTACATAATCCCCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTACACGTCAGCACGAAACT",
            b"GCTAAAGACAATTACATAATCCCCCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTACACGTCAGCACGAAACT",
        ];
        let truth: &[u8] =
            b"GCTAAAGACAATTACATAATCCCCCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGTTACACGTCAGCACGAAACT";
        let c = crate::consensus(&reads, 0, &PoaConfig::default()).unwrap();
        assert_eq!(
            c.sequence,
            truth,
            "got {}",
            String::from_utf8_lossy(&c.sequence)
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

#[cfg(test)]
mod tests_homopolymer_probe {
    use super::*;

    fn read_fasta(path: &str) -> Vec<Vec<u8>> {
        let fa = std::fs::read_to_string(path).expect("fasta");
        let mut reads = Vec::new();
        let mut cur = Vec::new();
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
        reads
    }

    /// Dump MSA columns of a real/handmade read set to see whether a single
    /// interrupting base between two variable homopolymers lands in one column
    /// or is split across several. PROBE_FA overrides the input.
    #[test]
    #[ignore = "probe; run with --ignored --nocapture"]
    fn probe_split_interrupt_column() {
        let path = std::env::var("PROBE_FA").unwrap_or_else(|_| {
            concat!(env!("CARGO_MANIFEST_DIR"), "/homopolymers/msa.fasta").into()
        });
        let reads = read_fasta(&path);
        let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let mut order: Vec<usize> = (0..reads.len()).collect();
        order.sort_by_key(|&i| reads[i].len());
        let med = order[order.len() / 2];
        let mut g = Poa::new(&reads[med], PoaConfig::default());
        for (i, r) in reads.iter().enumerate() {
            if i != med {
                g.add_read(r);
            }
        }
        let (node_col, ncol) = g.msa_columns().expect("acyclic columns");
        let mut cols: Vec<Vec<usize>> = vec![Vec::new(); ncol];
        for (nd, &c) in node_col.iter().enumerate() {
            cols[c].push(nd);
        }
        for (c, nodes) in cols.iter().enumerate() {
            let s: Vec<String> = nodes
                .iter()
                .map(|&n| {
                    format!(
                        "n{}:{} cov{} del{}",
                        n, g.nodes[n].base as char, g.nodes[n].cov, g.nodes[n].del
                    )
                })
                .collect();
            eprintln!("col {c:4}: {}", s.join(" | "));
        }
        eprintln!("heaviest : {}", String::from_utf8_lossy(&g.consensus()));
        eprintln!(
            "majority : {}",
            String::from_utf8_lossy(&g.consensus_majority())
        );
        eprintln!(
            "best-fit : {}",
            String::from_utf8_lossy(&g.consensus_best_fit(&refs))
        );
    }

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn unit(&mut self) -> f64 {
            self.next() as f64 / (1u64 << 31) as f64
        }
    }

    /// Run-length jitter: 0 with p0, else +/-1 (mostly) or +/-2.
    fn jitter(rng: &mut Lcg, p0: f64) -> i32 {
        let u = rng.unit();
        if u < p0 {
            return 0;
        }
        let mag = if rng.unit() < 0.8 { 1 } else { 2 };
        if rng.unit() < 0.5 { -mag } else { mag }
    }

    fn rand_flank(rng: &mut Lcg, n: usize) -> Vec<u8> {
        // No A at the junctions so the runs are unambiguous.
        let mut v: Vec<u8> = (0..n).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect();
        let k = v.len();
        v[k - 1] = b'C';
        v[0] = b'T';
        v
    }

    /// Sweep: flankL + A^a G A^b + flankR, each read jitters a and b
    /// independently. Count how often the interrupt G is split across MSA
    /// columns, and how often the public `consensus()` output is wrong
    /// (G lost/duplicated, or run lengths off).
    #[test]
    #[ignore = "probe; run with --ignored --nocapture"]
    fn probe_interrupt_sweep() {
        for &(a0, b0) in &[(12usize, 11usize), (8, 8), (20, 20)] {
            for &p0 in &[0.7, 0.5, 0.3] {
                for &depth in &[10usize, 20, 40] {
                    let trials = 200;
                    let (mut split, mut g_wrong, mut len_wrong) = (0, 0, 0);
                    for t in 0..trials {
                        let mut rng = Lcg(1 + t as u64 * 7919 + depth as u64 * 31);
                        let lf = rand_flank(&mut rng, 30);
                        let mut rf = rand_flank(&mut rng, 30);
                        rf[0] = b'C';
                        let mk = |a: usize, b: usize| {
                            let mut r = lf.clone();
                            r.extend(std::iter::repeat_n(b'A', a));
                            r.push(b'G');
                            r.extend(std::iter::repeat_n(b'A', b));
                            r.extend_from_slice(&rf);
                            r
                        };
                        let reads: Vec<Vec<u8>> = (0..depth)
                            .map(|_| {
                                let a = (a0 as i32 + jitter(&mut rng, p0)).max(1) as usize;
                                let b = (b0 as i32 + jitter(&mut rng, p0)).max(1) as usize;
                                mk(a, b)
                            })
                            .collect();
                        let truth = mk(a0, b0);
                        let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
                        // Column split check on the same graph the library builds.
                        let mut order: Vec<usize> = (0..reads.len()).collect();
                        order.sort_by_key(|&i| reads[i].len());
                        let med = order[order.len() / 2];
                        let mut g = Poa::new(&reads[med], PoaConfig::default());
                        for (i, r) in reads.iter().enumerate() {
                            if i != med {
                                g.add_read(r);
                            }
                        }
                        let (node_col, _) = g.msa_columns().expect("acyclic");
                        let mut gcols: Vec<usize> = (0..g.nodes.len())
                            .filter(|&n| g.nodes[n].base == b'G' && g.nodes[n].cov >= 2)
                            .map(|n| node_col[n])
                            .collect();
                        // Only the interrupt region: columns past the left flank.
                        gcols.retain(|&c| c >= lf.len() && c < lf.len() + a0 + b0 + 6);
                        gcols.sort_unstable();
                        gcols.dedup();
                        if gcols.len() > 1 {
                            split += 1;
                        }
                        let cons = crate::consensus(&refs, 0, &PoaConfig::default())
                            .unwrap()
                            .sequence;
                        let mid =
                            &cons[lf.len().min(cons.len())..cons.len().saturating_sub(rf.len())];
                        let ng = mid.iter().filter(|&&c| c == b'G').count();
                        if ng != 1 {
                            g_wrong += 1;
                        } else if cons != truth {
                            len_wrong += 1;
                        }
                    }
                    eprintln!(
                        "A{a0}GA{b0} p0={p0} depth={depth:2}: G split across columns {split:3}/{trials} | consensus G lost/dup {g_wrong:3} | runs wrong {len_wrong:3}"
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod tests_homopolymer_compare {
    use crate::config::PoaConfig;

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn unit(&mut self) -> f64 {
            self.next() as f64 / (1u64 << 31) as f64
        }
    }

    fn jitter(rng: &mut Lcg, p0: f64) -> i32 {
        if rng.unit() < p0 {
            return 0;
        }
        let mag = if rng.unit() < 0.8 { 1 } else { 2 };
        if rng.unit() < 0.5 { -mag } else { mag }
    }

    /// Raw vs homopolymer-compressed consensus on flankL + A^12 G A^11 + flankR,
    /// with per-read run-length jitter and random substitutions (rate `e`).
    /// Reports exact-match rate and G-lost/duplicated rate per mode.
    #[test]
    #[ignore = "probe; run with --ignored --nocapture"]
    fn probe_hpc_vs_raw() {
        let (a0, b0) = (12usize, 11usize);
        for &e in &[0.0, 0.01, 0.03] {
            for &p0 in &[0.7, 0.5, 0.3] {
                for &depth in &[10usize, 20, 40] {
                    let trials = 200;
                    let modes = [
                        (0usize, false),
                        (2, false),
                        (4, false),
                        (2, true),
                        (4, true),
                    ];
                    let mut exact = [0usize; 5];
                    let mut gbad = [0usize; 5];
                    for t in 0..trials {
                        let mut rng = Lcg(1 + t as u64 * 7919 + depth as u64 * 31);
                        let flank = |rng: &mut Lcg| -> Vec<u8> {
                            (0..30).map(|_| b"CGT"[(rng.next() % 3) as usize]).collect()
                        };
                        let lf = flank(&mut rng);
                        let rf = flank(&mut rng);
                        let mk = |a: usize, b: usize| {
                            let mut r = lf.clone();
                            r.extend(std::iter::repeat_n(b'A', a));
                            r.push(b'G');
                            r.extend(std::iter::repeat_n(b'A', b));
                            r.extend_from_slice(&rf);
                            r
                        };
                        let truth = mk(a0, b0);
                        let reads: Vec<Vec<u8>> = (0..depth)
                            .map(|_| {
                                let a = (a0 as i32 + jitter(&mut rng, p0)).max(1) as usize;
                                let b = (b0 as i32 + jitter(&mut rng, p0)).max(1) as usize;
                                let mut r = mk(a, b);
                                for x in r.iter_mut() {
                                    if rng.unit() < e {
                                        let alts: Vec<u8> =
                                            b"ACGT".iter().copied().filter(|&c| c != *x).collect();
                                        *x = alts[(rng.next() % 3) as usize];
                                    }
                                }
                                r
                            })
                            .collect();
                        let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
                        for (mi, &(m, best)) in modes.iter().enumerate() {
                            let mut cfg = PoaConfig {
                                homopolymer_min_run: m,
                                ..PoaConfig::default()
                            };
                            if std::env::var("MF").is_ok() {
                                cfg.consensus_mode =
                                    crate::config::ConsensusMode::MajorityFrequency;
                            }
                            let cons = if m == 0 || best {
                                crate::consensus(&refs, 0, &cfg).unwrap().sequence
                            } else {
                                crate::hpc::consensus_compressed(&refs, &cfg).sequence
                            };
                            if cons == truth {
                                exact[mi] += 1;
                            } else if std::env::var("SHOW").is_ok()
                                && m == 4
                                && best
                                && e > 0.02
                                && p0 > 0.6
                                && depth == 20
                                && t < 15
                            {
                                let f = |s: &[u8]| {
                                    String::from_utf8_lossy(&s[25..s.len().saturating_sub(25)])
                                        .into_owned()
                                };
                                eprintln!("  truth {}\n  k2    {}", f(&truth), f(&cons));
                            }
                            // The interrupt region is the only place with A's (flanks are C/G/T).
                            let first = cons.iter().position(|&c| c == b'A');
                            let last = cons.iter().rposition(|&c| c == b'A');
                            let ng = match (first, last) {
                                (Some(f), Some(l)) => {
                                    cons[f..=l].iter().filter(|&&c| c == b'G').count()
                                }
                                _ => 0,
                            };
                            if ng != 1 {
                                gbad[mi] += 1;
                            }
                        }
                    }
                    let fmt = |v: &[usize; 5]| {
                        v.iter()
                            .map(|x| format!("{x:3}"))
                            .collect::<Vec<_>>()
                            .join(" ")
                    };
                    eprintln!(
                        "e={e:.2} p0={p0} d={depth:2} | exact raw/k2/k4/k2best/k4best: {} | G bad: {}",
                        fmt(&exact),
                        fmt(&gbad)
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod tests_hp_methods {
    use crate::config::PoaConfig;

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn unit(&mut self) -> f64 {
            self.next() as f64 / (1u64 << 31) as f64
        }
    }

    fn jitter(rng: &mut Lcg, p0: f64) -> i32 {
        if rng.unit() < p0 {
            return 0;
        }
        let mag = if rng.unit() < 0.8 { 1 } else { 2 };
        if rng.unit() < 0.5 { -mag } else { mag }
    }

    fn flank(rng: &mut Lcg) -> Vec<u8> {
        let mut v: Vec<u8> = (0..30).map(|_| b"CGT"[(rng.next() % 3) as usize]).collect();
        v[0] = b'C';
        v[29] = b'C';
        v
    }

    fn methods() -> Vec<(&'static str, PoaConfig)> {
        let d = PoaConfig::default();
        vec![
            ("raw", d.clone()),
            // "!" = engine with gap costs directly, no splice/fit guard.
            (
                "gap!",
                PoaConfig {
                    homopolymer_gap: Some((-2, -1)),
                    ..d.clone()
                },
            ),
            (
                "gap",
                PoaConfig {
                    homopolymer_gap: Some((-2, -1)),
                    ..d.clone()
                },
            ),
            (
                "hpc2",
                PoaConfig {
                    homopolymer_min_run: 2,
                    ..d.clone()
                },
            ),
            (
                "hpc4",
                PoaConfig {
                    homopolymer_min_run: 4,
                    ..d.clone()
                },
            ),
            (
                "hpc2+gap",
                PoaConfig {
                    homopolymer_min_run: 2,
                    homopolymer_gap: Some((-2, -1)),
                    ..d.clone()
                },
            ),
        ]
    }

    /// `make(rng) -> (truth, reads)`; prints exact-match count and mean edit
    /// distance per method.
    /// (truth, reads)
    type Case = (Vec<u8>, Vec<Vec<u8>>);

    fn run(label: &str, trials: usize, make: &dyn Fn(&mut Lcg) -> Case) {
        let ms = methods();
        let mut exact = vec![0usize; ms.len()];
        let mut edits = vec![0usize; ms.len()];
        let mut units = vec![0usize; ms.len()];
        for t in 0..trials {
            let mut rng = Lcg(17 + t as u64 * 104729);
            let (truth, reads) = make(&mut rng);
            let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
            for (mi, (_, cfg)) in ms.iter().enumerate() {
                let cons = if ms[mi].0.ends_with('!') {
                    crate::engine_consensus(&refs, cfg).sequence
                } else {
                    crate::consensus(&refs, 0, cfg).unwrap().sequence
                };
                let ed = crate::poa2::tests::lev(&cons, &truth);
                edits[mi] += ed;
                let nt = |s: &[u8]| s.iter().filter(|&&c| c == b'T').count();
                if nt(&cons) == nt(&truth) {
                    units[mi] += 1;
                }
                if ed == 0 {
                    exact[mi] += 1;
                }
            }
        }
        let cells: Vec<String> = ms
            .iter()
            .enumerate()
            .map(|(i, (n, _))| {
                format!(
                    "{n} {:3}/{:3} ({:.2})",
                    exact[i],
                    units[i],
                    edits[i] as f64 / trials as f64
                )
            })
            .collect();
        eprintln!("{label} | {}", cells.join(" | "));
    }

    fn subs(rng: &mut Lcg, r: &mut [u8], e: f64) {
        for x in r.iter_mut() {
            if rng.unit() < e {
                let alts: Vec<u8> = b"ACGT".iter().copied().filter(|&c| c != *x).collect();
                *x = alts[(rng.next() % 3) as usize];
            }
        }
    }

    /// Per method: exact matches / T-count matches (unit count for the AT
    /// repeat) out of N, and mean edit distance.
    #[test]
    #[ignore = "probe; run with --ignored --nocapture"]
    fn probe_hp_methods() {
        let trials: usize = std::env::var("TRIALS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(100);
        let only = std::env::var("ONLY").unwrap_or_default();
        for &e in &[0.0, 0.01, 0.03] {
            for &p0 in &[0.7, 0.4] {
                for &depth in &[10usize, 20, 40] {
                    if only.is_empty() || only == "g" {
                        run(
                            &format!("G   e={e:.2} p0={p0} d={depth:2}"),
                            trials,
                            &|rng| {
                                let (lf, rf) = (flank(rng), flank(rng));
                                let mk = |a: usize, b: usize| {
                                    let mut r = lf.clone();
                                    r.extend(std::iter::repeat_n(b'A', a));
                                    r.push(b'G');
                                    r.extend(std::iter::repeat_n(b'A', b));
                                    r.extend_from_slice(&rf);
                                    r
                                };
                                let reads = (0..depth)
                                    .map(|_| {
                                        let a = (12 + jitter(rng, p0)).max(1) as usize;
                                        let b = (11 + jitter(rng, p0)).max(1) as usize;
                                        let mut r = mk(a, b);
                                        subs(rng, &mut r, e);
                                        r
                                    })
                                    .collect();
                                (mk(12, 11), reads)
                            },
                        );
                    }
                    if only.is_empty() || only == "r" {
                        run(
                            &format!("AT  e={e:.2} p0={p0} d={depth:2}"),
                            trials,
                            &|rng| {
                                let (lf, rf) = (flank(rng), flank(rng));
                                // 15 units of A{n}T; ~30% of units are genuinely A3/A5
                                // in every read (shared), the rest A4.
                                let units: Vec<usize> = (0..15)
                                    .map(|_| {
                                        let u = rng.unit();
                                        if u < 0.15 {
                                            3
                                        } else if u < 0.3 {
                                            5
                                        } else {
                                            4
                                        }
                                    })
                                    .collect();
                                let mk = |ns: &[usize]| {
                                    let mut r = lf.clone();
                                    for &n in ns {
                                        r.extend(std::iter::repeat_n(b'A', n));
                                        r.push(b'T');
                                    }
                                    r.extend_from_slice(&rf);
                                    r
                                };
                                let reads = (0..depth)
                                    .map(|_| {
                                        // Per-read unique noise on every A run.
                                        let ns: Vec<usize> = units
                                            .iter()
                                            .map(|&n| (n as i32 + jitter(rng, p0)).max(1) as usize)
                                            .collect();
                                        let mut r = mk(&ns);
                                        subs(rng, &mut r, e);
                                        r
                                    })
                                    .collect();
                                (mk(&units), reads)
                            },
                        );
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests_support_quality {
    use super::*;

    #[test]
    fn support_phred_values() {
        // 10 reads agree, no alternative: P = 2^-11.
        assert_eq!(support_phred(10, 0), 33);
        // 6-4 split: P(Bin(11, 1/2) > 6) = 562 / 2048.
        assert_eq!(support_phred(6, 4), 6);
        // tie, and no spanning reads at all: P = 1/2.
        assert_eq!(support_phred(5, 5), 3);
        assert_eq!(support_phred(0, 0), 3);
        // the alternative outvotes the call
        assert_eq!(support_phred(2, 8), 0);
        // deep agreement is capped
        assert_eq!(support_phred(200, 0), MAX_SUPPORT_QUALITY);
        assert_eq!(support_phred(10_000, 3), MAX_SUPPORT_QUALITY);
    }

    #[test]
    fn support_phred_is_monotone() {
        for a in 0..40 {
            for m in 0..10 {
                assert!(support_phred(a + 1, m) >= support_phred(a, m));
                assert!(support_phred(a, m + 1) <= support_phred(a, m));
            }
        }
    }
}

#[cfg(test)]
mod tests_left_aligned_gaps {
    //! The traceback runs right to left and prefers to keep matching, so a gap
    //! opens only when forced: at the LEFTMOST of its equivalent positions
    //! inside a run or tandem repeat. Every read therefore places the same
    //! event at the same graph position, which is what lines reads up in the
    //! graph (and what an explicit left-normalisation pass would do; one was
    //! tried and moved nothing on the plain engine). These tests guard that
    //! property against changes to the traceback order.
    use super::*;

    /// True if no indel run in `ops` can slide one base left: the op before a
    /// run is not an exact match whose base equals the run's last base.
    fn left_aligned(g: &Poa, read: &[u8], ops: &[Op]) -> bool {
        let mut rp = 0usize;
        let mut prev_match: Option<(usize, u8)> = None; // (node, read base)
        let mut i = 0;
        while i < ops.len() {
            match ops[i] {
                Op::Match(t) => {
                    prev_match = Some((t, read[rp]));
                    rp += 1;
                    i += 1;
                }
                Op::Ins(_) | Op::Del(_) => {
                    let is_del = matches!(ops[i], Op::Del(_));
                    let mut e = i;
                    while e < ops.len() && matches!(ops[e], Op::Del(_)) == is_del {
                        if !is_del && !matches!(ops[e], Op::Ins(_)) {
                            break;
                        }
                        if !is_del {
                            rp += 1;
                        }
                        e += 1;
                    }
                    if let Some((p, rb)) = prev_match {
                        let pb = g.nodes[p].base;
                        let last = match ops[e - 1] {
                            Op::Del(d) => g.nodes[d].base,
                            Op::Ins(b) => b,
                            Op::Match(_) => unreachable!(),
                        };
                        if pb == rb && last == pb {
                            return false;
                        }
                    }
                    prev_match = None;
                    i = e;
                }
            }
        }
        true
    }

    #[test]
    fn checker_rejects_a_right_placed_gap() {
        use Op::{Del as D, Match as M};
        let g = Poa::new(b"GCAAAAT", PoaConfig::default());
        // GCAAAT with the deleted A on the last A of the run, then on the first
        assert!(!left_aligned(
            &g,
            b"GCAAAT",
            &[M(0), M(1), M(2), M(3), M(4), D(5), M(6)]
        ));
        assert!(left_aligned(
            &g,
            b"GCAAAT",
            &[M(0), M(1), D(2), M(3), M(4), M(5), M(6)]
        ));
    }

    #[test]
    fn gaps_in_runs_and_repeats_are_left_aligned() {
        let cfg = PoaConfig {
            homopolymer_min_run: 0,
            homopolymer_gap: None,
            ..PoaConfig::default()
        };
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut rnd = |n: usize| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state % n as u64) as usize
        };
        let flank_l = b"GATCGTACGTTAGCATGCTG";
        let flank_r = b"TTGCAGTCAGGATCCATGCA";
        for trial in 0..200 {
            // a homopolymer or dinucleotide tract, then one indel inside it
            let tract: Vec<u8> = if trial % 2 == 0 {
                vec![b"ACGT"[rnd(4)]; 6 + rnd(10)]
            } else {
                b"CA".repeat(5 + rnd(10))
            };
            let mut seed = flank_l.to_vec();
            seed.extend(&tract);
            seed.extend(flank_r);
            let unit = if trial % 2 == 0 { 1 } else { 2 };
            let at = flank_l.len() + unit * rnd(tract.len() / unit - 1);
            let mut read = seed.clone();
            if rnd(2) == 0 {
                read.drain(at..at + unit);
            } else {
                let ins = read[at..at + unit].to_vec();
                read.splice(at..at, ins);
            }
            let g = Poa::new(&seed, cfg.clone());
            let (topo, rank) = g.topo_order();
            let ops = g.align(&read, &topo, &rank);
            assert!(
                left_aligned(&g, &read, &ops),
                "trial {trial}: gap not left-aligned\nseed {}\nread {}\nops {ops:?}",
                String::from_utf8_lossy(&seed),
                String::from_utf8_lossy(&read)
            );
        }
    }
}
