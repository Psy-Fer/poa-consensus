use crate::config::{AlignmentMode, ConsensusMode, PoaConfig};
use crate::error::PoaError;
use crate::types::{BubbleSite, Consensus, CoverageGap, GapKind, GraphStats};
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

/// Half-width of the per-node j-window given to a locked winning arm node.
/// Wide enough to absorb ~5 cumulative indels at 5% ONT error rate.
const LOCK_EPS: usize = 5;

/// Half-width of the local-population window used by `coverage_threshold`'s
/// absolute-floor fallback (see `local_population_profile` and its call site
/// in `consensus()`).
///
/// Empirically tuned, not derived from first principles -- both a much
/// smaller and a much larger radius were tried and rejected:
///
/// - A read-length-scaled radius (the intuitive first choice: a read can
///   only contribute evidence within roughly one read-length of where it
///   starts) recovers almost none of the true partial-population fix --
///   see the call site's comment.
/// - A *small* radius (in this file's usual small-scale-jitter range, e.g.
///   `LOCK_EPS = 5` / `MINI_EPS_BASE = 3`) gives the best recovery on the
///   partial-read scenario this fix targets, but regresses three
///   pre-existing majority-length-wins tests on real repeat data
///   (`diag_sca8_real_sequences_no_flank`,
///   `diag_dab1_sca37_attttc_lookahead_arm_length_bias`,
///   `structural_phasing_no_contamination_on_noisy_periodic_repeat`): in a
///   tandem repeat, a genuine minority's extra repeat units sit only a few
///   bases past the majority's true boundary, so a small window bridges
///   from the minority's own low coverage straight into the adjacent
///   majority peak and wrongly rescues it -- reproducing the exact bug
///   pattern (Known Bugs #1-#10) this crate has spent the most effort
///   fixing, just via a new mechanism.
/// - A sweep against those three regressions found their exact crossover
///   radius at 18-20 (`diag_sca8`) and 40-45 (the periodic-phasing test);
///   `LOCAL_POP_RADIUS = 50` sits safely above both with margin, confirmed
///   against the full test suite (zero regressions) rather than living
///   right at the measured edge.
/// - This is a real, measured trade-off, not a free lunch: at radius=50 the
///   partial-read regression scenario (`gen_short_reads_long_region_ont_r10`
///   in `bench/compare_callers.py --general`) improves from edit=352 to
///   edit=282 against truth (vs external tools' 250/234) -- better, but it
///   does not fully close the gap the way a smaller, unsafe radius would.
///   Closing it further needs a mechanism that distinguishes repeat-driven
///   column proximity from genuine partial-population sparsity, which a
///   single fixed window cannot do; out of scope here.
const LOCAL_POP_RADIUS: usize = 50;

// ─── Minimizer anchoring constants ───────────────────────────────────────────

/// k-mer length for spine/read minimizer seeding.
const MINI_K: usize = 15;
/// Minimizer window width: select the minimum-hash k-mer from each window of
/// this many consecutive k-mers.  Density ≈ 1/MINI_W anchors per base.
const MINI_W: usize = 10;
/// Minimum per-anchor j-window half-width (absorbs indel errors near anchors).
/// Covers ~3 cumulative indels; keeps the correct j inside the intersection when
/// the spine and read have small local alignment offsets at the anchor position.
const MINI_EPS_BASE: usize = 3;
/// Minimum anchor chain length required before applying any anchor refinement.
/// Chosen to be above the expected number of error-induced wrong anchors in a
/// repetitive read (~1-3 for typical ONT depths), but below the expected number
/// of correct anchors in a non-repetitive read (≈read_len × 0.46 / MINI_W for
/// 5% error / k=15; ~23 for a 500bp read).  Flank-only reads (200bp of flank
/// from 100bp ANCHOR_PAD × 2) yield ≈9 flank anchors — below this threshold,
/// so anchors are safely disabled when the spine is repetitive and only flank
/// k-mers would match.
const MINI_MIN_CHAIN: usize = 15;
// Minimum chain density (anchors per base × MINI_W) to trust the chain.
// In non-repetitive sequence, density ≈ 1.0 (one anchor per MINI_W bases).
// In repetitive sequence, error-induced coincidences produce a sparse chain
// (density << 0.5) that is dominated by wrong anchors and must be ignored.
// Threshold: chain_len * MINI_W * 2 >= read_len  (≡ density ≥ 0.5 / MINI_W).
// Equivalently: require chain_len >= read_len / (2 * MINI_W).

// ─── Align scratch ───────────────────────────────────────────────────────────

/// Per-call scratch reused across `align()` calls to avoid per-call heap
/// allocation for the lock-window tables.  Both vecs are sorted by node_idx
/// and binary-searched; they are empty when no lock fired, so the check is
/// O(1) in the common case.
struct AlignScratch {
    /// Locked winning arm nodes: (node_idx, j_center as u32).
    lock_node_j: Vec<(usize, u32)>,
    /// Locked exit nodes: (node_idx, j_center as u32).
    lock_exit_j: Vec<(usize, u32)>,
}

impl AlignScratch {
    fn new() -> Self {
        Self {
            lock_node_j: Vec::new(),
            lock_exit_j: Vec::new(),
        }
    }

    fn clear(&mut self) {
        self.lock_node_j.clear();
        self.lock_exit_j.clear();
    }

    fn get_node_j(&self, node_idx: usize) -> Option<usize> {
        self.lock_node_j
            .binary_search_by_key(&node_idx, |&(idx, _)| idx)
            .ok()
            .map(|pos| self.lock_node_j[pos].1 as usize)
    }

    fn get_exit_j(&self, node_idx: usize) -> Option<usize> {
        self.lock_exit_j
            .binary_search_by_key(&node_idx, |&(idx, _)| idx)
            .ok()
            .map(|pos| self.lock_exit_j[pos].1 as usize)
    }
}

/// Current weight of the edge `from -> to`, or 0 if no such edge exists yet.
///
/// Used to break DP score ties in favour of the already-better-supported
/// transition.  Without this, ties are broken by `in_edges` iteration order
/// (effectively read-arrival order), which lets content-identical repeat
/// registers (e.g. two phase-shifted encodings of the same homopolymer run)
/// split support roughly evenly across reads instead of converging onto one
/// path — silently fragmenting the graph in periodic sequence.
///
/// Returns the *total* (matched + deleted) traversal count, not matched-only.
/// This is a live-alignment tie-break ("which predecessor is more
/// established"), a different question from `heaviest_path`'s final
/// path-selection scoring. Open Question #2 (design/graph_data_model_rework.md
/// Phase 1): tried both matched-only and total here, A/B against the full
/// `cargo test` suite and `bench/validate.py` -- neither discriminated between
/// them (every currently-passing case stayed passing, and the 3 scenarios
/// `bench/validate.py` fails are pre-existing, present identically at the
/// unmodified `HEAD` in an isolated worktree, unrelated to this choice). With
/// no evidence favouring a change, kept at total weight -- unchanged from
/// before the Match/Delete edge-weight split.
#[inline]
fn edge_weight(nodes: &[Node], from: usize, to: usize) -> i32 {
    nodes[from]
        .out_edges
        .iter()
        .find(|&&(t, _)| t == to)
        .map(|&(_, ew)| ew.total())
        .unwrap_or(0)
}

#[inline]
fn safe_add(a: i32, b: i32) -> i32 {
    if a == UNSET {
        UNSET
    } else {
        a.saturating_add(b)
    }
}

// ─── Graph types ─────────────────────────────────────────────────────────────

/// Per-edge traversal weight, split by how the read traversed it.
///
/// A single conflated `i32` weight cannot tell "N reads confirmed this base"
/// (`Match`) apart from "N reads skipped past whatever's here" (`Delete`) --
/// and both increment the *same* edge when they share a predecessor. See
/// `design/graph_data_model_rework.md` for the audit that motivated this
/// split (confirmed on real RFC1 data: nodes reached by one pure-delete
/// in-edge and one pure-match in-edge to the same target).
///
/// `Insert` traffic is folded into `matched`, not tracked separately: once a
/// node exists (whether created by `Insert` or already present), every later
/// read that agrees with its base routes through `Match`, so `Insert`'s own
/// founding traversal is accounting-identical to a `Match` -- both mean "a
/// read confirms this exact base occurs here." A three-way split would add a
/// bucket that never needs to be read differently from `matched`.
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
struct EdgeWeight {
    /// Match traversals, plus the founding Insert that created the target node.
    matched: i32,
    /// Delete traversals (read skipped past the target node without consuming a base).
    deleted: i32,
}

impl EdgeWeight {
    fn total(&self) -> i32 {
        self.matched + self.deleted
    }
}

struct Node {
    base: u8,
    /// (successor_node_index, edge_weight)
    out_edges: Vec<(usize, EdgeWeight)>,
    /// predecessor node indices
    in_edges: Vec<usize>,
    /// reads that produced a Match op at this node (not Delete ops)
    coverage: u32,
    /// reads that produced a Delete op at this node (traversed without consuming a base)
    delete_count: u32,
    /// Cached `(fork, arm_entry)`: the nearest ancestor with 2+ out-edges (a
    /// real fork) and the specific child of that fork this node's own
    /// single-successor chain descends from -- `None` if no such ancestor
    /// exists (this node has never had a fork anywhere back along its
    /// creation chain). `arm_entry` is cached alongside `fork` because the
    /// interior filter's plurality comparison needs to know *which* of the
    /// fork's out-edges this node's arm belongs to, not just that a fork
    /// exists.
    ///
    /// Populated incrementally in `add_to_graph`/`propagate_fork_if_new` --
    /// see design/graph_data_model_rework.md Phase 4. Can go *stale in one
    /// direction only*: pointing farther back than the true nearest fork,
    /// never invalid (a node that is a fork never stops being one, since
    /// out-edges are only ever added, and `arm_entry` is only ever set
    /// alongside a `fork` that was genuinely a fork when set). Reconvergence
    /// points (2+ in-edges) and anything past them are deliberately left
    /// untouched by propagation -- ambiguous which incoming arm "owns" them,
    /// so whatever was established at their own creation stays authoritative.
    nearest_fork: Option<(usize, usize)>,
}

pub struct PoaGraph {
    nodes: Vec<Node>,
    config: PoaConfig,
    /// number of reads integrated (seed + subsequent)
    n_reads: usize,
    /// original reads stored for per-allele sub-graph reconstruction
    reads: Vec<Vec<u8>>,
    /// (from_node, to_node) → sorted list of read indices that genuinely
    /// confirmed this edge's target (Match, or the founding Insert that
    /// created it). Does NOT include reads that merely deleted through the
    /// target. Mirrors the `coverage`/`delete_count` idiom at the node level
    /// (see design/graph_data_model_rework.md Phase 2).
    edge_reads: HashMap<(usize, usize), Vec<u32>>,
    /// **Vestigial under pure bypass (revised Phase 1 of
    /// design/bypass_edge_delete_rework.md).** Formerly `(from, to)` → read
    /// indices that Deleted through an edge; it was always write-only (no
    /// production consumer ever read it -- confirmed by exhaustive grep). Under
    /// pure bypass a deleting read touches the skipped node only via
    /// `delete_count`, never via an edge, so this map is no longer populated at
    /// all; it stays initialised-empty. Left in place rather than removed to
    /// keep this phase's delta minimal; removable in a later cleanup.
    #[allow(dead_code)]
    edge_delete_reads: HashMap<(usize, usize), Vec<u32>>,
    /// Per-node bypass edges: from-node index → list of `(to_node, weight)`.
    /// A bypass edge records that a read skipped one or more nodes via a run
    /// of consecutive `Delete` ops and then reconnected -- by a clean Match --
    /// onto a node already in the graph (the run's entry predecessor bypasses
    /// straight to that resume node). This is a *separate* structural edge
    /// around the skipped node(s) -- exactly how abPOA and SPOA represent a
    /// deletion (see `design/architecture_comparison_abpoa_spoa.md` Finding 1,
    /// and `design/bypass_edge_delete_rework.md` for the full phased plan).
    ///
    /// **Pure-bypass representation (revised Phase 1).** The skipped nodes are
    /// touched only by `delete_count`; no laundered matched through-edge is
    /// created (that is the delta from the superseded dual-bookkeeping Phase 1
    /// commit). A resume that instead *diverges into new structure*
    /// (mismatch/insert) is recorded as an ordinary minority `out_edge`, not a
    /// bypass -- so bypass edges here are exactly "rejoined the existing main
    /// path after a skip."
    ///
    /// **Phase 1 status: written but not yet read by any consumer.**
    /// `heaviest_path` (the intended reader) does not consult it until Phase 2;
    /// the interior filter, `find_bubbles`/`find_structural_bubbles`,
    /// `compute_stats`, and `verify_reuse_chain` never will (they iterate
    /// `Node.out_edges`, which a bypass edge is deliberately kept out of).
    ///
    /// Weight is a plain `i32`, deliberately *not* an [`EdgeWeight`]: a bypass
    /// edge is definitionally delete-traffic, and keeping it in this separate
    /// structure (never in `out_edges`) is what makes it invisible to every
    /// `matched`-based out-edge consumer for free -- so it cannot manufacture a
    /// false bubble, at any deletion rate.
    bypass_edges: HashMap<usize, Vec<(usize, i32)>>,
    /// number of times the long-unbanded warning was emitted
    warnings: usize,
    /// Cached heaviest-path spine, recomputed adaptively rather than every read.
    cached_spine: Vec<(usize, u8, i32)>,
    /// `n_reads` at the time the spine was last recomputed.
    spine_updated_at: usize,
    /// Current recompute interval; doubles when the spine is stable, resets
    /// to 1 when it changes significantly.
    spine_interval: usize,
    /// Reusable scratch for lock-window tables; avoids per-call heap allocation.
    align_scratch: AlignScratch,
    /// Minimizer index over the cached spine sequence; rebuilt whenever
    /// `cached_spine` is refreshed.  Maps k-mer hash → spine rank (index into
    /// `cached_spine`).  Only hashes that appear exactly once in the spine are
    /// stored — duplicate k-mers cannot serve as unambiguous anchors.
    spine_mers: HashMap<u64, u32>,
    /// Per-fork content-addressable arm index: fork node index → (full
    /// characterized edit's bases → existing arm's start node index).
    /// Populated in `add_to_graph` whenever a new divergence arm is created;
    /// consulted before creating a new arm so an exact-duplicate edit (same
    /// fork, same complete base sequence) reuses the existing node chain
    /// instead of fragmenting into a new one. Scoped per-fork (not global) so
    /// reuse can only ever happen among arms of the *same* divergence point;
    /// see design/graph_data_model_rework.md Phase 3. Does NOT solve fuzzy
    /// near-duplicates (edits with the same effect but different incidental
    /// length/shape) -- only byte-identical edits hash-match.
    fork_arm_index: HashMap<usize, HashMap<Vec<u8>, usize>>,
    /// Set when any `add_read` call needed `align_with_retry`'s pass 2 or 3
    /// (i.e. the configured band was too narrow for at least one read).
    /// Read by `build_graph` (`src/lib.rs`): a graph built from a *mix* of
    /// some reads succeeding on the configured (narrow) band and others
    /// only succeeding after a wider retry is not the same as building the
    /// whole graph with a single, consistent band from the start -- in a
    /// periodic/repetitive locus, different reads can settle on different,
    /// individually-plausible diagonals (Known Bug #3's mechanism),
    /// fragmenting bubble structure that a uniformly-unbanded build would
    /// not. `build_graph` uses this flag to rebuild the whole graph
    /// unbanded from scratch when set, rather than trusting a
    /// mixed-band graph as-is.
    used_band_retry: bool,
}

// ─── DP cell ─────────────────────────────────────────────────────────────────

/// Sentinel topo index meaning "came from the virtual start node".
const VIRTUAL: u32 = u32::MAX;

#[derive(Clone, Copy, PartialEq, Eq)]
enum State {
    M,
    I,
    D,
}

/// DP cell: 8 bytes (i32 score + u32 pred_t, no padding).
/// pred_t stores the topo-order row of the best predecessor, or VIRTUAL.
#[derive(Clone, Copy)]
struct Cell {
    score: i32,
    pred_t: u32,
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
        nearest_fork: None,
    });
    idx
}

/// Bounded depth for forward propagation when a predecessor transitions from
/// a single out-edge to a fork (2+ out-edges), fixing up the *pre-existing*
/// arm's descendants whose cached `nearest_fork` now points farther back than
/// this newly-forked predecessor. Mirrors `ARM_MAX_DEPTH`'s existing bound
/// (the same order of magnitude already accepted elsewhere in this file for
/// walking a single-successor chain) rather than inventing a new, unproven
/// limit -- see design/graph_data_model_rework.md Phase 4 for why an
/// *unbounded* walk here would repeat the exact class of mistake Phase 3
/// found (nodes/walks with no depth limit blowing up on long unbranched runs).
const FORK_PROPAGATE_MAX_DEPTH: usize = ARM_MAX_DEPTH;

/// Called immediately after a new out-edge is added from `from` to some
/// node, whether that node is brand new or already existed. If this is
/// `from`'s *second* out-edge (a fresh 1 -> 2 transition, i.e. `from` just
/// became a fork), walks forward along `from`'s *original* out-edge (the one
/// that existed before this new one -- `out_edges[0]`) updating every
/// descendant's cached `nearest_fork` to `Some((from, orig_target))`, until
/// hitting a node that is itself already a fork (anything past it already
/// has a closer fork of its own to be attributed to) or a reconvergence
/// point (2+ in-edges -- ambiguous which incoming arm "owns" it, left
/// untouched), or the depth bound. A no-op if `from` was already a fork
/// before this edge (its pre-existing arms' descendants are already
/// correctly attributed) or is still not one (this was its first edge).
fn propagate_fork_if_new(nodes: &mut [Node], from: usize) {
    if nodes[from].out_edges.len() != 2 {
        return;
    }
    let orig_target = nodes[from].out_edges[0].0;
    let mut cur = orig_target;
    for _ in 0..FORK_PROPAGATE_MAX_DEPTH {
        if nodes[cur].in_edges.len() > 1 {
            break;
        }
        nodes[cur].nearest_fork = Some((from, orig_target));
        match nodes[cur].out_edges.as_slice() {
            [(next, _)] => cur = *next,
            _ => break,
        }
    }
}

/// Sets a freshly-created node's own `nearest_fork`. Does NOT propagate to
/// `p`'s *pre-existing* arm -- that has to be deferred until the current
/// read's entire traceback has been processed (see `add_to_graph`'s
/// `to_propagate` list), not run eagerly here. Confirmed empirically (not
/// assumed) that eager, mid-read propagation is wrong: a single read's own
/// later ops can turn a node the propagation already updated into a
/// reconvergence point (2+ in-edges) a few steps later in that *same*
/// read's own processing -- at the moment of the eager update, that
/// downstream node still looks like an ordinary single-predecessor node, so
/// the "don't touch reconvergence points" guard in `propagate_fork_if_new`
/// can't see the problem coming. Running propagation only after the whole
/// read is processed avoids this: by then, its own reconvergence, if any, has
/// already happened, so `in_edges.len()` reflects the read's final effect.
fn set_new_node_own_fork(nodes: &mut [Node], p: usize, new_idx: usize) {
    if nodes[p].out_edges.len() >= 2 {
        // p is a fork (whether it just became one or already was) --
        // new_idx starts a fresh arm directly off it.
        nodes[new_idx].nearest_fork = Some((p, new_idx));
    } else {
        nodes[new_idx].nearest_fork = nodes[p].nearest_fork;
    }
}

/// Adds a brand-new edge with weight 1, attributed as genuine confirmation
/// (`matched: 1, deleted: 0`). Used only for the seed's initial linear chain,
/// where the seed read's own bases are by definition confirmed. Calls
/// `propagate_fork_if_new` directly (not deferred): the seed's entire chain
/// is built in one uninterrupted pass with exactly one out-edge per node, so
/// there is no "this read's own later ops" to race against -- the deferred-
/// propagation subtlety only applies within a single `add_to_graph` call.
fn add_edge(nodes: &mut [Node], from: usize, to: usize) {
    nodes[from].out_edges.push((
        to,
        EdgeWeight {
            matched: 1,
            deleted: 0,
        },
    ));
    nodes[to].in_edges.push(from);
    propagate_fork_if_new(nodes, from);
}

/// Increments an existing edge's weight, or creates it at weight 1, routing
/// the `+1` to `.matched` or `.deleted` depending on how this read traversed
/// it (`is_delete`). Returns `true` if a genuinely new edge was created (as
/// opposed to an existing one being incremented) -- callers in `add_to_graph`
/// use this to defer `propagate_fork_if_new(nodes, from)` until the current
/// read's traceback has been fully processed, not run it immediately.
#[must_use]
fn increment_or_add_edge(nodes: &mut [Node], from: usize, to: usize, is_delete: bool) -> bool {
    debug_assert_ne!(
        from, to,
        "increment_or_add_edge: refusing to create/increment a literal self-loop"
    );
    for (succ, ew) in nodes[from].out_edges.iter_mut() {
        if *succ == to {
            if is_delete {
                ew.deleted += 1;
            } else {
                ew.matched += 1;
            }
            return false;
        }
    }
    let ew = if is_delete {
        EdgeWeight {
            matched: 0,
            deleted: 1,
        }
    } else {
        EdgeWeight {
            matched: 1,
            deleted: 0,
        }
    };
    nodes[from].out_edges.push((to, ew));
    nodes[to].in_edges.push(from);
    true
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

/// Number of `node`'s in-edges whose own source-side weight (matched-only, the
/// same axis `find_structural_bubbles`/`find_bubbles` already score arms on)
/// clears `sig_threshold`. Used by [`materialize_arm_len_tolerant`] to tell a
/// genuine reconvergence (2+ real paths meeting) apart from one real path plus
/// incidental noise (e.g. one read's stray Delete-turned-edge, or a single
/// misaligned base creating a low-weight side edge).
fn real_in_edge_count(nodes: &[Node], node: usize, sig_threshold: i32) -> usize {
    nodes[node]
        .in_edges
        .iter()
        .filter(|&&p| {
            nodes[p]
                .out_edges
                .iter()
                .find(|&&(to, _)| to == node)
                .is_some_and(|&(_, ew)| ew.matched >= sig_threshold)
        })
        .count()
}

/// Like [`materialize_arm_len`], but tolerant of minor internal noise: a fork
/// or reconvergence only ends the arm when it is a *real* one -- i.e. when
/// 2+ of the edges involved individually clear `sig_threshold` -- rather than
/// on any branch at all, no matter how small.
///
/// # Why this exists (see `find_structural_bubbles`)
///
/// A genuinely long structural difference between two alleles (e.g. GAA×30
/// vs GAA×100) does not, in real data, materialise as one perfectly clean
/// unbranched chain of matches for its entire length: ordinary per-read
/// sequencing noise creates small side-branches and reconvergences scattered
/// throughout it (a single read's substitution error creating a 1-weight
/// alternate node; a single read's Delete creating a low-weight direct edge
/// past a node others still Match through). The strict, non-tolerant walk
/// stops at the *first* such disruption, so it measures only the distance to
/// the first read's noise, not the arm's true span -- confirmed empirically
/// on real data (RFC1 hap2's whole "delete-driven forkless gap" class of
/// known bugs is this same phenomenon in a different codepath) and on a
/// synthetic 150 bp non-repetitive structural-insertion control at realistic
/// (~6%) per-read error rates, where the strict walk also undercounts span
/// for a perfectly clean, unambiguous, non-periodic difference.
///
/// # What counts as "minor" and why the bound is safe
///
/// A branch or reconvergence is "minor" exactly when it does NOT clear
/// `sig_threshold` -- the identical `min_allele_freq`-derived vote count
/// `find_bubbles`/`find_structural_bubbles` already use to decide "is this
/// arm significant enough to be a competing-allele candidate" everywhere
/// else in this file. This is deliberate, not a new, separately-tuned bar:
/// anything this function is willing to walk through is, by the SAME
/// definition already governing every other bubble decision in this module,
/// not significant enough to be treated as a second real path. It cannot
/// silently swallow a genuine second haplotype nested inside the arm being
/// measured, because a genuine one clears the same bar used to find it in
/// the first place -- if 2+ of the candidate edges at a fork (or 2+ of a
/// node's in-edges) each clear `sig_threshold` independently, that is by
/// construction a real branch/reconvergence, and the walk stops there exactly
/// as the strict version would.
///
/// Bounded to `max_depth` steps, identical to the caller's existing bound
/// (`phasing_bubble_min_span + 1`): the caller only ever needs to know
/// whether the span reaches `phasing_bubble_min_span`, so once the walk has
/// gone that far without hitting a REAL fork or reconvergence the answer is
/// already known, and there is no reason -- or additional risk -- in
/// continuing further. This keeps the same bounded-walk shape as this file's
/// existing patterns (`ARM_MAX_DEPTH`, `FORK_SEARCH_HOPS`) rather than
/// introducing an open-ended search; unlike reusing the whole-graph
/// `compute_bubble_ranges` path-budget algorithm (which tracks true
/// reconvergence of *every* out-edge, including noise ones, and can mark an
/// entire open-ended tail of the graph as "one giant unresolved bubble" when
/// a noise branch never reconverges before the graph ends), this walk only
/// ever looks at the single arm being measured and never grows unbounded.
fn materialize_arm_len_tolerant(
    nodes: &[Node],
    start: usize,
    max_depth: usize,
    sig_threshold: i32,
) -> usize {
    let mut len = 0;
    let mut cur = start;
    for _ in 0..max_depth {
        // After the first step, only a REAL reconvergence (2+ in-edges each
        // individually clearing sig_threshold) ends the arm.
        if len >= 1 && real_in_edge_count(nodes, cur, sig_threshold) > 1 {
            break;
        }
        len += 1;
        match nodes[cur].out_edges.as_slice() {
            [] => break,
            [(next, _)] => cur = *next,
            edges => {
                let real_arms = edges
                    .iter()
                    .filter(|&&(_, ew)| ew.matched >= sig_threshold)
                    .count();
                if real_arms > 1 {
                    break; // a genuine nested fork -- this arm ends here.
                }
                // At most one edge is significant; the rest are noise. A
                // significant edge's weight is by definition >= sig_threshold
                // and thus >= any noise edge's weight, so max_by_key always
                // picks it when one exists, and picks an arbitrary (but
                // deterministic) noise edge when none do.
                let (next, _) = edges.iter().max_by_key(|&&(_, ew)| ew.matched).unwrap();
                cur = *next;
            }
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
/// For each topo rank, compute the bubble it belongs to as `Some((entry_t, exit_t))`,
/// or `None` for spine nodes outside any bubble.
///
/// Uses a generation-counter scratch (`gen_mark`) so the inner bubble search never
/// allocates: each new bubble increments the generation and marks nodes with it
/// instead of resetting a bool array.
fn compute_bubble_ranges(nodes: &[Node], topo: &[usize]) -> Vec<Option<(usize, usize)>> {
    let n = topo.len();
    let nn = nodes.len();
    let mut ranges: Vec<Option<(usize, usize)>> = vec![None; n];
    // gen_mark[node_idx] == current_gen  ⟺  node is inside the active bubble search.
    // Allocated once; never cleared between bubbles.
    let mut gen_mark = vec![0u32; nn];
    let mut current_gen = 0u32;

    let mut t = 0;
    while t < n {
        if nodes[topo[t]].out_edges.len() >= 2 {
            // Start a new bubble search from entry topo rank t.
            current_gen = current_gen.wrapping_add(1);
            let entry_node = topo[t];
            gen_mark[entry_node] = current_gen;
            let mut outstanding = nodes[entry_node].out_edges.len();
            let mut exit_t = None;

            for (tt, &node_idx) in topo.iter().enumerate().skip(t + 1) {
                let bubble_in = nodes[node_idx]
                    .in_edges
                    .iter()
                    .filter(|&&p| gen_mark[p] == current_gen)
                    .count();
                if bubble_in == 0 {
                    continue;
                }
                gen_mark[node_idx] = current_gen;
                outstanding -= bubble_in;
                if outstanding == 0 {
                    exit_t = Some(tt);
                    break;
                }
                outstanding += nodes[node_idx].out_edges.len();
            }

            if let Some(et) = exit_t {
                ranges[t..=et].fill(Some((t, et)));
                t = et + 1;
            } else {
                // Open-ended bubble: mark everything remaining conservatively.
                ranges[t..n].fill(Some((t, n - 1)));
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
/// 4096 covers AAGGG × 800 (one arm) and similar long STR bubbles.
const ARM_MAX_DEPTH: usize = 4096;

/// Walks one bubble arm from `start_idx` forward, collecting node indices
/// until the bubble exit (topo rank `exit_t`) is reached.  The exit node
/// itself is NOT included — it runs its own windowed DP normally.
///
/// Returns `None` if the arm has internal branching (a nested bubble) or
/// exceeds the depth cap; both are signals to fall back to windowed DP.
/// Walk one bubble arm from `start_idx` to (but not including) the exit node,
/// appending node indices into `out`.  Clears `out` first.
/// Returns `false` if the arm has an internal branch or exceeds `ARM_MAX_DEPTH`.
fn collect_arm_nodes(
    nodes: &[Node],
    rank_of: &[usize],
    start_idx: usize,
    exit_t: usize,
    out: &mut Vec<usize>,
) -> bool {
    out.clear();
    let mut cur = start_idx;
    for _ in 0..ARM_MAX_DEPTH {
        if rank_of[cur] == exit_t {
            return true;
        }
        out.push(cur);
        match nodes[cur].out_edges.as_slice() {
            [(next, _)] => cur = *next,
            _ => return false,
        }
    }
    false
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

// ─── Slide-and-lock ──────────────────────────────────────────────────────────

/// Number of consecutive positions where one arm uniquely matches the query
/// before we commit to it.  Two positions guards against single-base read
/// errors flipping the decision in repetitive sequence.
const SLIDE_MIN_CONSEC: usize = 2;

/// Slide all arms against `query[j_entry..]` simultaneously, base by base,
/// until one arm can be uniquely identified as the correct path.
///
/// **Lock conditions** (first to fire wins):
/// - *Exhaustion*: an arm runs out of nodes while the query has remaining
///   bases and at least one other arm can still continue.  The exhausted arm
///   cannot account for the remaining query and is eliminated.
/// - *Unique mismatch*: exactly one alive arm matches the query base at the
///   current position for `SLIDE_MIN_CONSEC` consecutive steps.  All others
///   are eliminated.
/// - *Single survivor*: only one arm remains after eliminations above.
///
/// Returns the winning arm index into `all_arms`, or `None` if no lock fires
/// within the arm bounds.  `None` falls back to windowed DP, which is always
/// correct.
///
/// **No retroactive gap-state fill is required.**  Because slide-and-lock is
/// decided at the bubble *entry* node (before arm nodes are processed in topo
/// order), the winning arm's nodes proceed through the normal windowed DP and
/// accumulate correct M/I/D cells for traceback.
fn slide_lock(
    all_arms: &[Vec<usize>],
    nodes: &[Node],
    query: &[u8],
    j_entry: usize,
) -> Option<usize> {
    let n_arms = all_arms.len();
    let remaining = query.len().saturating_sub(j_entry);
    if remaining == 0 || n_arms < 2 {
        return None;
    }

    // Only slide when ALL arms are at least LOOKAHEAD_K nodes long.
    // Short arms are produced by individual read errors (1-4 nodes):
    //   - insertion bubbles have one empty arm (0 nodes, direct edge to exit)
    //     which exhausts immediately, incorrectly eliminating the no-insertion path
    //   - 2-3 consecutive substitution errors produce 2-3 node arms where a
    //     single read error can produce SLIDE_MIN_CONSEC false unique matches
    // Requiring all arms ≥ LOOKAHEAD_K ensures we only slide on structural
    // bubbles of at least one full repeat unit.
    let min_arm_len = all_arms.iter().map(|a| a.len()).min().unwrap_or(0);
    if min_arm_len < LOOKAHEAD_K {
        return None;
    }

    let max_steps = all_arms
        .iter()
        .map(|a| a.len())
        .max()
        .unwrap_or(0)
        .min(remaining);

    let mut alive = vec![true; n_arms];
    let mut alive_count = n_arms;
    let mut consec_unique = 0usize;
    let mut consec_winner: Option<usize> = None;

    for i in 0..max_steps {
        // ── Exhaustion ────────────────────────────────────────────────────────
        // Eliminate arms that have run out of nodes while at least one other arm
        // can still continue.  If ALL alive arms exhaust simultaneously the query
        // continues past the bubble at the shared exit — no lock from exhaustion.
        let any_continuing = (0..n_arms).any(|idx| alive[idx] && all_arms[idx].len() > i);
        if any_continuing {
            for idx in 0..n_arms {
                if alive[idx] && all_arms[idx].len() <= i {
                    alive[idx] = false;
                    alive_count -= 1;
                }
            }
        }
        if alive_count <= 1 {
            break;
        }

        // ── Unique mismatch ───────────────────────────────────────────────────
        let q = query[j_entry + i];
        // Bounds check on arm index: alive arms always have len > i here (ensured
        // by the exhaustion block above which eliminates len <= i arms first).
        let matched: Vec<bool> = (0..n_arms)
            .map(|idx| alive[idx] && all_arms[idx].len() > i && nodes[all_arms[idx][i]].base == q)
            .collect();
        let match_count = matched.iter().filter(|&&m| m).count();
        let mismatch_count = (0..n_arms)
            .filter(|&idx| alive[idx] && !matched[idx])
            .count();

        if match_count == 1 && mismatch_count > 0 {
            let candidate = matched.iter().position(|&m| m).unwrap();
            if consec_winner == Some(candidate) {
                consec_unique += 1;
            } else {
                consec_unique = 1;
                consec_winner = Some(candidate);
            }
            if consec_unique >= SLIDE_MIN_CONSEC {
                for idx in 0..n_arms {
                    if alive[idx] && !matched[idx] {
                        alive[idx] = false;
                        alive_count -= 1;
                    }
                }
                break;
            }
        } else {
            consec_unique = 0;
            consec_winner = None;
        }
    }

    if alive_count == 1 {
        alive.iter().position(|&a| a)
    } else {
        None
    }
}

// ─── Minimizer anchoring ─────────────────────────────────────────────────────

#[inline]
fn encode_base(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Compute minimizers for `seq`: for each sliding window of `w` consecutive
/// k-mers, select the (hash, start-position) pair with the minimum hash.
/// Returns deduplicated (hash, seq_position) pairs in order of appearance.
/// K-mers containing non-ACGT bases are skipped (treated as invalid).
fn compute_minimizers(seq: &[u8], k: usize, w: usize) -> Vec<(u64, usize)> {
    let n = seq.len();
    if n < k {
        return vec![];
    }
    let n_kmers = n - k + 1;
    let mask = if k * 2 < 64 {
        (1u64 << (k * 2)) - 1
    } else {
        u64::MAX
    };

    // Rolling 2-bit kmer hash; None at positions containing non-ACGT bases.
    let mut kmer_hashes: Vec<Option<u64>> = Vec::with_capacity(n_kmers);
    let mut hash: u64 = 0;
    let mut valid: usize = 0;
    for (i, &b) in seq.iter().enumerate() {
        match encode_base(b) {
            Some(bits) => {
                hash = ((hash << 2) | bits) & mask;
                valid += 1;
            }
            None => {
                hash = 0;
                valid = 0;
            }
        }
        if i + 1 >= k {
            kmer_hashes.push(if valid >= k { Some(hash) } else { None });
        }
    }

    // Slide a window of w k-mer positions; pick the minimum-hash position.
    // win = w.min(n_kmers) ensures the loop is safe when the sequence is short.
    let win = w.min(n_kmers);
    let mut result: Vec<(u64, usize)> = Vec::new();
    let mut last: Option<(u64, usize)> = None;
    for start in 0..=(n_kmers - win) {
        let mut best_hash = u64::MAX;
        let mut best_pos = start;
        for (off, entry) in kmer_hashes[start..start + win].iter().enumerate() {
            if let Some(h) = entry {
                if *h < best_hash {
                    best_hash = *h;
                    best_pos = start + off;
                }
            }
        }
        if best_hash == u64::MAX {
            continue;
        }
        let entry = (best_hash, best_pos);
        if last != Some(entry) {
            last = Some(entry);
            result.push(entry);
        }
    }
    result
}

/// Build the minimizer index for the cached spine.
/// Maps k-mer hash → spine rank (index in `spine`).
/// Only k-mers that appear exactly once in the spine are kept.
fn build_spine_mers(spine: &[(usize, u8, i32)], k: usize, w: usize) -> HashMap<u64, u32> {
    let spine_seq: Vec<u8> = spine.iter().map(|&(_, b, _)| b).collect();
    let mers = compute_minimizers(&spine_seq, k, w);
    // Count occurrences and record the first-seen spine rank per hash.
    let mut counts: HashMap<u64, (u32, u32)> = HashMap::new();
    for (hash, pos) in mers {
        let e = counts.entry(hash).or_insert((0, pos as u32));
        e.0 += 1;
    }
    counts
        .into_iter()
        .filter(|(_, (c, _))| *c == 1)
        .map(|(h, (_, p))| (h, p))
        .collect()
}

/// Return the indices into `vals` that form the longest strictly increasing
/// subsequence.  O(n log n) patience-sort with predecessor backtracking.
fn lis_indices(vals: &[usize]) -> Vec<usize> {
    let n = vals.len();
    if n == 0 {
        return vec![];
    }
    let mut tails: Vec<usize> = Vec::new();
    let mut tail_idx: Vec<usize> = Vec::new();
    let mut pred: Vec<Option<usize>> = vec![None; n];
    for i in 0..n {
        let v = vals[i];
        let p = tails.partition_point(|&t| t < v);
        if p == tails.len() {
            tails.push(v);
            tail_idx.push(i);
        } else {
            tails[p] = v;
            tail_idx[p] = i;
        }
        pred[i] = if p > 0 { Some(tail_idx[p - 1]) } else { None };
    }
    let lis_len = tails.len();
    let mut chain = Vec::with_capacity(lis_len);
    let mut cur = Some(tail_idx[lis_len - 1]);
    while let Some(i) = cur {
        chain.push(i);
        cur = pred[i];
    }
    chain.reverse();
    chain
}

/// Build the colinear anchor chain for one read against the current spine.
///
/// Matches each read minimizer against `spine_mers` to produce candidate
/// `(read_pos, topo_rank)` pairs, then runs LIS to keep the longest
/// monotonically-increasing sub-sequence (guaranteeing colinearity).
/// Returns the chain sorted by topo_rank.
fn build_anchor_chain(
    read_mers: &[(u64, usize)],
    spine_mers: &HashMap<u64, u32>,
    spine: &[(usize, u8, i32)],
    rank_of: &[usize],
) -> Vec<(usize, usize)> {
    let mut candidates: Vec<(usize, usize)> = read_mers
        .iter()
        .filter_map(|&(hash, read_pos)| {
            spine_mers.get(&hash).map(|&sr| {
                let node_idx = spine[sr as usize].0;
                (read_pos, rank_of[node_idx])
            })
        })
        .collect();
    if candidates.is_empty() {
        return vec![];
    }
    // Sort by topo_rank; LIS on read_pos then gives a colinear monotone chain.
    candidates.sort_unstable_by_key(|&(_, t)| t);
    let read_positions: Vec<usize> = candidates.iter().map(|&(r, _)| r).collect();
    let idxs = lis_indices(&read_positions);
    idxs.iter().map(|&i| candidates[i]).collect()
}

/// Compute the anchor-derived j-window for topo row `t`.
///
/// Returns `Some((lo, hi))` only when `t` is the exact topo rank of an anchor.
/// Interpolation between anchors is intentionally omitted: in repetitive
/// sequence, error k-mers at different repeat positions can produce coincident
/// k-mer hashes (same error type, same cycle position → identical k-mer), which
/// would create wrong anchors that are colinear with the correct flank anchors
/// and survive the LIS filter.  An exact-match-only policy is safe because a
/// wrong anchor centered far from the correct j produces an empty intersection
/// with the existing window → the existing window is used unchanged.
fn anchor_j_bounds(
    t: usize,
    anchors: &[(usize, usize)], // (read_pos, topo_rank), sorted by topo_rank
    l: usize,
) -> Option<(usize, usize)> {
    if anchors.is_empty() {
        return None;
    }
    // Binary search for an exact topo_rank match.
    let pos = anchors.partition_point(|&(_, at)| at <= t);
    if pos > 0 && anchors[pos - 1].1 == t {
        let r = anchors[pos - 1].0;
        let lo = r.saturating_sub(MINI_EPS_BASE).max(1);
        let hi = (r + MINI_EPS_BASE).min(l);
        Some((lo, hi))
    } else {
        None
    }
}

/// Intersect `(j_lo, j_hi)` with the anchor-derived window for topo row `t`.
///
/// Only applied to spine nodes that are outside bubble regions and not already
/// constrained by slide-lock windows.  Returns the original window unchanged
/// when anchors provide no constraint or the intersection is empty (wrong
/// anchor — correctness fallback to existing window).
///
/// Also skips the anchor when its upper bound falls below `j_center`
/// (j_pred_max + 1): a shifted anchor that would exclude the expected
/// diagonal causes the diagonal skip to write outside the valid band,
/// silently discarding the correct alignment path.
#[allow(clippy::too_many_arguments)]
#[inline]
fn anchor_refine_spine(
    j_lo: usize,
    j_hi: usize,
    t: usize,
    anchors: &[(usize, usize)],
    l: usize,
    is_source: bool,
    is_spine: bool,
    in_bubble: bool,
    is_locked: bool,
    j_center: usize,
) -> (usize, usize) {
    if anchors.is_empty() || is_source || !is_spine || in_bubble || is_locked {
        return (j_lo, j_hi);
    }
    match anchor_j_bounds(t, anchors, l) {
        Some((alo, ahi)) => {
            // A shifted anchor whose upper bound sits below the expected diagonal
            // would exclude the correct alignment path.  Discard it.
            if ahi < j_center {
                return (j_lo, j_hi);
            }
            let nlo = j_lo.max(alo);
            let nhi = j_hi.min(ahi);
            if nlo <= nhi { (nlo, nhi) } else { (j_lo, j_hi) }
        }
        None => (j_lo, j_hi),
    }
}

// ─── Bubble-aware DP alignment ────────────────────────────────────────────────

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
// Bounds-safe read from a banded matrix (m or ins).  j=0 is always UNSET
// for these two matrices.  Returns UNSET when j falls outside [j_lo, j_hi].
#[inline(always)]
fn gs(mat: &[Cell], t: usize, j: usize, j_lo: usize, j_hi: usize, rw: usize) -> i32 {
    if j == 0 || j < j_lo || j > j_hi {
        UNSET
    } else {
        mat[t * rw + (j - j_lo)].score
    }
}

// Bounds-safe read for the del matrix.  j=0 is stored in del0 (separate array).
#[inline(always)]
fn gsd(
    del: &[Cell],
    del0: &[Cell],
    t: usize,
    j: usize,
    j_lo: usize,
    j_hi: usize,
    rw: usize,
) -> i32 {
    if j == 0 {
        del0[t].score
    } else if j < j_lo || j > j_hi {
        UNSET
    } else {
        del[t * rw + (j - j_lo)].score
    }
}

#[allow(clippy::too_many_arguments)]
fn align(
    nodes: &[Node],
    topo: &[usize],
    rank_of: &[usize],
    spine: &[(usize, u8, i32)],
    query: &[u8],
    cfg: &PoaConfig,
    scratch: &mut AlignScratch,
    anchors: &[(usize, usize)], // (read_pos, topo_rank) colinear chain from minimizer seeding
) -> Result<Vec<AlignOp>, PoaError> {
    scratch.clear();
    let n = topo.len();
    let l = query.len();

    // band_width=0 with adaptive_band=false is genuinely unbanded (O(V*L) DP,
    // see spine_margin below) -- warn once per call on long reads so callers
    // aren't surprised by the memory/time cost.  Never blocks the call; set
    // PoaConfig::warn_on_long_unbanded = false to suppress.
    if cfg.warn_on_long_unbanded && cfg.band_width == 0 && !cfg.adaptive_band && l > 1000 {
        eprintln!(
            "poa-consensus: warning: unbanded alignment (band_width=0, \
             adaptive_band=false) on a {l} bp read -- this scales as \
             O(read_len * graph_len); consider a banded or adaptive PoaConfig \
             for large graphs, or set warn_on_long_unbanded=false to suppress"
        );
    }

    let go = cfg.gap_open;
    let ge = cfg.gap_extend;
    let semi = cfg.alignment_mode == AlignmentMode::SemiGlobal;

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

    // Spine margin: matches the DP band width so memory usage scales with the
    // configured band rather than always allocating SPINE_MARGIN_MIN columns.
    // Uses the same adaptive formula as main (b + f*L).  band_width=0 with
    // adaptive_band=false is the documented "unbanded / full NW over DAG" case
    // (see PoaConfig::band_width): spine_margin is set to the full query length
    // so row_width below covers [1, l] for every row, with no artificial cap.
    // This is deliberately expensive at long read lengths (see CLAUDE.md scale
    // table) -- callers who explicitly ask for unbanded get real O(V*L) DP, not
    // a silently narrow fallback band.
    let spine_margin: usize = if cfg.adaptive_band {
        let w = cfg.adaptive_band_b + (cfg.adaptive_band_f * l as f32).ceil() as usize;
        let w = if cfg.band_width > 0 {
            w.max(cfg.band_width)
        } else {
            w
        };
        w.max(4)
    } else if cfg.band_width > 0 {
        cfg.band_width
    } else {
        l
    };

    // Banded DP tables.  Each row stores j = j_lo_arr[t]..=j_hi_arr[t] (j >= 1).
    // Spine rows use ±spine_margin around the diagonal; bubble rows widen by the
    // bubble span.  The j=0 del column is stored separately in del0 because
    // m[t][0] and ins[t][0] are structurally UNSET and never stored.
    // +2 ensures the left-edge j-1 access from the next row stays in range.
    let row_width = (2 * spine_margin + 2).min(l).max(1);

    let mut del0 = vec![Cell::unset(); n];
    let bsize = n * row_width;
    let mut m = vec![Cell::unset(); bsize];
    let mut ins = vec![Cell::unset(); bsize];
    let mut del = vec![Cell::unset(); bsize];
    // Per-row band: j_lo_arr[t] is the lowest stored j (>= 1).
    // j_hi_arr[t] is set when the row is processed; reads before then return UNSET.
    let mut j_lo_arr = vec![1usize; n];
    let mut j_hi_arr = vec![0usize; n];

    // Reusable scratch for bubble-arm materialisation (avoids per-bubble alloc).
    let mut arm_scratch: Vec<usize> = Vec::new();
    let mut all_arms: Vec<Vec<usize>> = Vec::new();

    for (t, &node_idx) in topo.iter().enumerate() {
        // Lookahead: skip nodes on losing bubble arms entirely.
        if lookahead_skip[node_idx] {
            continue;
        }

        let node_base = nodes[node_idx].base;
        let is_source = nodes[node_idx].in_edges.is_empty();

        // ── j range: compute window FIRST (needed for diagonal skip write) ─────
        // Spine node  : j_center = predecessor's best_j + 1  (±SPINE_MARGIN)
        // Bubble node : j_center = bubble-entry's predecessor j, hi widened by span
        // Source node : clamped full range [1, row_width] (large-j cells are UNSET
        //               in predecessors so they're never winners; safe to clamp).
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
            (1usize, row_width.min(l))
        } else if let Some(j_center) = scratch.get_node_j(node_idx) {
            // Locked winning arm node: tight per-node window centred at the
            // exact query column this depth corresponds to.
            let lo = j_center.saturating_sub(LOCK_EPS).max(1);
            let hi = (j_center + LOCK_EPS).min(l);
            (lo, hi)
        } else if let Some(j_center) = scratch.get_exit_j(node_idx) {
            // Locked exit node: spine-width window centred at winning arm's
            // terminal query column, so it reads the arm's actual best cells.
            let lo = j_center.saturating_sub(spine_margin).max(1);
            let hi = (j_center + spine_margin).min(l);
            (lo, hi)
        } else {
            match bubble_ranges[t] {
                Some((entry_t, exit_t)) => {
                    if t == entry_t {
                        bubble_entry_j[entry_t] = j_pred_max;
                    }
                    let bej = bubble_entry_j[entry_t];
                    let bubble_span = exit_t.saturating_sub(entry_t) + 1;
                    let lo = bej.saturating_sub(spine_margin).max(1);
                    let hi = (bej + bubble_span.min(spine_margin) + spine_margin).min(l);
                    (lo, hi)
                }
                None => {
                    let j_center = j_pred_max.saturating_add(1);
                    let lo = j_center.saturating_sub(spine_margin).max(1);
                    let hi = (j_center + spine_margin).min(l);
                    (lo, hi)
                }
            }
        };
        // Clamp hi to the banded row width.
        let j_hi = j_hi.min(j_lo + row_width - 1);
        // Narrow with minimizer anchor constraint for spine nodes that are
        // outside bubble regions and not already pinned by slide-lock windows.
        // Never widens the window; falls back to (j_lo, j_hi) when the anchor
        // constraint is absent, would produce an empty interval, or would
        // exclude the expected diagonal (shifted-anchor safety guard).
        let (j_lo, j_hi) = anchor_refine_spine(
            j_lo,
            j_hi,
            t,
            anchors,
            l,
            is_source,
            on_spine[node_idx],
            bubble_ranges[t].is_some(),
            scratch.get_node_j(node_idx).is_some() || scratch.get_exit_j(node_idx).is_some(),
            j_pred_max.saturating_add(1),
        );
        j_lo_arr[t] = j_lo;
        j_hi_arr[t] = j_hi;

        // ── Diagonal skip ──────────────────────────────────────────────────────
        // O(1) fast path for consecutive exact matches along the spine.
        //
        // Extended to fire at bubble entries: if a spine node has multiple
        // out-edges (e.g. a low-weight error arm), a 1-base pre-resolve checks
        // whether query[bj+1] uniquely matches the spine successor's first base.
        // If so, losing arms are marked lookahead_skip and the skip fires as if
        // the node had a single out-edge.  Error arms never block spine sliding.
        //
        // In-edge guard is relaxed: predecessors already in lookahead_skip (from
        // a previously resolved bubble) are ignored, so exit nodes of resolved
        // bubbles re-enter the skip path immediately.
        #[cfg(test)]
        if !is_source {
            NODE_COUNTER.with(|c| c.set(c.get() + 1));
        }

        if !is_source && on_spine[node_idx] {
            if let Some(prev_sp) = spine_prev[node_idx] {
                // Relaxed in-edge check: fast for the common single-predecessor case;
                // falls through to a loop only when there are extra predecessors that
                // might all be lookahead_skip (exit nodes of resolved bubbles).
                let active_pred_ok = match nodes[node_idx].in_edges.as_slice() {
                    [p] => *p == prev_sp,
                    _ => {
                        nodes[node_idx]
                            .in_edges
                            .iter()
                            .all(|&p| p == prev_sp || lookahead_skip[p])
                            && nodes[node_idx].in_edges.contains(&prev_sp)
                    }
                };

                if active_pred_ok {
                    let pt = rank_of[prev_sp];
                    let bj = best_j_per_t[pt];
                    // bj+1 must land inside this row's own band, or the write at
                    // line ~1472 (`bj + 1 - j_lo`) underflows: j_lo here is this
                    // node's window, which can start after bj+1 (e.g. a locked
                    // per-node window centred elsewhere), independent of the
                    // predecessor's diagonal position.
                    if bj < l && bj + 1 >= j_lo && bj < j_hi && node_base == query[bj] {
                        let m_prev = gs(&m, pt, bj, j_lo_arr[pt], j_hi_arr[pt], row_width);
                        let i_prev = gs(&ins, pt, bj, j_lo_arr[pt], j_hi_arr[pt], row_width);
                        let d_prev =
                            gsd(&del, &del0, pt, bj, j_lo_arr[pt], j_hi_arr[pt], row_width);
                        let pred_is_source = nodes[topo[pt]].in_edges.is_empty();
                        let i_ok = i_prev == UNSET || (pred_is_source && m_prev > i_prev);
                        let d_ok = d_prev == UNSET || (pred_is_source && m_prev > d_prev);
                        if m_prev != UNSET && i_ok && d_ok {
                            // For a bubble entry: 1-base pre-resolve against query[bj+1].
                            // Mark losing arms by setting only their first node as skip;
                            // the UNSET cascade through the arm is correct and avoids the
                            // alloc inside collect_arm_nodes on the hot path.
                            let do_skip = match nodes[node_idx].out_edges.as_slice() {
                                [_] => true,
                                _ if bj + 1 < l => {
                                    let next_q = query[bj + 1];
                                    let mut spine_succ = None;
                                    let mut resolved = true;
                                    for &(s, _) in &nodes[node_idx].out_edges {
                                        if on_spine[s] && spine_prev[s] == Some(node_idx) {
                                            if nodes[s].base == next_q {
                                                spine_succ = Some(s);
                                            } else {
                                                resolved = false;
                                                break;
                                            }
                                        } else if !lookahead_skip[s] && nodes[s].base == next_q {
                                            resolved = false;
                                            break;
                                        }
                                    }
                                    if resolved {
                                        if let Some(ss) = spine_succ {
                                            for &(s, _) in &nodes[node_idx].out_edges {
                                                if s != ss && !lookahead_skip[s] {
                                                    lookahead_skip[s] = true;
                                                    // Walk the rest of this losing arm, not
                                                    // just its first node: a multi-node
                                                    // insertion (e.g. a 2+ base indel) left
                                                    // every node past the first one neither
                                                    // on-spine nor marked, so it fell through
                                                    // to real windowed DP on every subsequent
                                                    // read, and -- worse -- kept the arm's
                                                    // reconvergence node's active_pred_ok
                                                    // check failing (it's still a live,
                                                    // unresolved in-edge), forcing the entire
                                                    // rest of that read into full DP too.
                                                    // Stops at the spine (correct
                                                    // reconvergence) or a further branch
                                                    // (left to its own resolution).
                                                    let mut cur = s;
                                                    for _ in 0..ARM_MAX_DEPTH {
                                                        match nodes[cur].out_edges.as_slice() {
                                                            [(next, _)]
                                                                if !on_spine[*next]
                                                                    && !lookahead_skip[*next] =>
                                                            {
                                                                lookahead_skip[*next] = true;
                                                                cur = *next;
                                                            }
                                                            _ => break,
                                                        }
                                                    }
                                                }
                                            }
                                            true
                                        } else {
                                            false
                                        }
                                    } else {
                                        false
                                    }
                                }
                                _ => false,
                            };

                            if do_skip {
                                let score = m_prev + cfg.match_score;
                                // bj+1 is the centre of this row's band: within [j_lo, j_hi].
                                m[t * row_width + (bj + 1 - j_lo)] = Cell {
                                    score,
                                    pred_t: pt as u32,
                                };
                                best_j_per_t[t] = bj + 1;
                                #[cfg(test)]
                                SKIP_COUNTER.with(|c| c.set(c.get() + 1));
                                continue;
                            }
                        }
                    }
                }
            }
        }

        // ── j = 0: delete-only ──────────────────────────────────────────────────
        // m[t][0] and ins[t][0] are structurally UNSET (never stored); only del[t][0]
        // (del0[t]) needs initialisation here.
        {
            let (mut best, mut best_pred) = (UNSET, 0u32);
            if is_source {
                let val = go + ge;
                if val > best {
                    best = val;
                    best_pred = VIRTUAL;
                }
            }
            for &p in &nodes[node_idx].in_edges {
                let pt = rank_of[p];
                // m[pt][0] is always UNSET; only del0 carries the j=0 chain.
                let vd = safe_add(del0[pt].score, ge);
                if vd > best {
                    best = vd;
                    best_pred = pt as u32;
                }
            }
            if best != UNSET {
                del0[t] = Cell {
                    score: best,
                    pred_t: best_pred,
                };
            }
        }

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
            let ixcur = t * row_width + (j - j_lo);

            // M[t][j]
            {
                let (mut best, mut best_pred) = (UNSET, VIRTUAL);
                // Free-start at j=1: applies when j_lo==1 (source or bubble nodes).
                if j == 1 && (is_source || semi) && sc > best {
                    best = sc;
                    best_pred = VIRTUAL;
                }
                if is_source && j > 1 {
                    let val = safe_add(go + (j as i32 - 1) * ge, sc);
                    if val != UNSET && val > best {
                        best = val;
                        best_pred = VIRTUAL;
                    }
                }
                // Ties across predecessors prefer the heavier-weighted edge (see
                // `edge_weight`); ties within the same predecessor keep the
                // Match > Insert > Delete order via the pre-existing `>` checks.
                // MAX so a pre-loop free-start `best` is never displaced by a tie.
                let mut best_edge_w = i32::MAX;
                for &p in &nodes[node_idx].in_edges {
                    let pt = rank_of[p];
                    let ew = edge_weight(nodes, p, node_idx);
                    let vm = safe_add(gs(&m, pt, j - 1, j_lo_arr[pt], j_hi_arr[pt], row_width), sc);
                    if vm != UNSET && (vm > best || (vm == best && ew > best_edge_w)) {
                        best = vm;
                        best_pred = pt as u32;
                        best_edge_w = ew;
                    }
                    let vi = safe_add(
                        gs(&ins, pt, j - 1, j_lo_arr[pt], j_hi_arr[pt], row_width),
                        sc,
                    );
                    if vi != UNSET && (vi > best || (vi == best && ew > best_edge_w)) {
                        best = vi;
                        best_pred = pt as u32;
                        best_edge_w = ew;
                    }
                    let vd = safe_add(
                        gsd(
                            &del,
                            &del0,
                            pt,
                            j - 1,
                            j_lo_arr[pt],
                            j_hi_arr[pt],
                            row_width,
                        ),
                        sc,
                    );
                    if vd != UNSET && (vd > best || (vd == best && ew > best_edge_w)) {
                        best = vd;
                        best_pred = pt as u32;
                        best_edge_w = ew;
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

            // I[t][j] — reads from same row at j-1.
            // When j == j_lo, j-1 is outside the banded window → treat as UNSET.
            {
                let (mut best, mut best_pred) = (UNSET, VIRTUAL);
                if is_source && j == 1 {
                    let val = go + ge;
                    if val > best {
                        best = val;
                        best_pred = VIRTUAL;
                    }
                }
                if j > j_lo {
                    let ixprev = t * row_width + (j - 1 - j_lo);
                    let vm = safe_add(m[ixprev].score, go + ge);
                    if vm != UNSET && vm > best {
                        best = vm;
                        best_pred = t as u32;
                    }
                    let vi = safe_add(ins[ixprev].score, ge);
                    if vi != UNSET && vi > best {
                        best = vi;
                        best_pred = t as u32;
                    }
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
                let (mut best, mut best_pred) = (UNSET, 0u32);
                // See the M[t][j] block above: ties across predecessors prefer
                // the heavier-weighted edge instead of in_edges iteration order.
                let mut best_edge_w = i32::MAX;
                for &p in &nodes[node_idx].in_edges {
                    let pt = rank_of[p];
                    let ew = edge_weight(nodes, p, node_idx);
                    let vm = safe_add(
                        gs(&m, pt, j, j_lo_arr[pt], j_hi_arr[pt], row_width),
                        go + ge,
                    );
                    if vm != UNSET && (vm > best || (vm == best && ew > best_edge_w)) {
                        best = vm;
                        best_pred = pt as u32;
                        best_edge_w = ew;
                    }
                    let vi = safe_add(
                        gs(&ins, pt, j, j_lo_arr[pt], j_hi_arr[pt], row_width),
                        go + ge,
                    );
                    if vi != UNSET && (vi > best || (vi == best && ew > best_edge_w)) {
                        best = vi;
                        best_pred = pt as u32;
                        best_edge_w = ew;
                    }
                    let vd = safe_add(
                        gsd(&del, &del0, pt, j, j_lo_arr[pt], j_hi_arr[pt], row_width),
                        ge,
                    );
                    if vd != UNSET && (vd > best || (vd == best && ew > best_edge_w)) {
                        best = vd;
                        best_pred = pt as u32;
                        best_edge_w = ew;
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

        // ── Lookahead resolve + slide-and-lock ───────────────────────────────
        // At a bubble entry node, attempt to identify the correct arm before
        // processing its nodes.  Losing arm nodes are marked in `lookahead_skip`
        // and their windowed DP is entirely skipped; the winning arm proceeds
        // through normal windowed DP, producing correct M/I/D cells for
        // traceback without any retroactive gap-state fill.
        //
        // Dispatch order:
        //   1. Lookahead  — score first LOOKAHEAD_K bases; fast, reliable for
        //                   structural bubbles (interruption motifs, allele SNPs)
        //                   where one arm is clearly better within one repeat unit.
        //   2. Slide-and-lock — slide base-by-base past LOOKAHEAD_K until an arm
        //                   exhausts or uniquely matches for SLIDE_MIN_CONSEC steps;
        //                   handles long repetitive bubbles (RFC1, C9orf72) where
        //                   arms share a long identical prefix.
        //   3. Windowed DP fallback — always correct; used when arms are complex
        //                   (internally branched) or no lock fires.
        if let Some((entry_t, exit_t)) = bubble_ranges[t] {
            if t == entry_t && nodes[node_idx].out_edges.len() >= 2 {
                let j_entry = best_j_per_t[t];
                if j_entry < l {
                    // Materialise each arm.  Arms with internal branches (nested
                    // bubbles) or exceeding ARM_MAX_DEPTH return None → complex=true
                    // → fall through to windowed DP.
                    all_arms.clear();
                    let mut complex = false;
                    for &(arm_start, _) in &nodes[node_idx].out_edges {
                        if collect_arm_nodes(nodes, rank_of, arm_start, exit_t, &mut arm_scratch) {
                            all_arms.push(arm_scratch.clone());
                        } else {
                            complex = true;
                            break;
                        }
                    }

                    if !complex && all_arms.len() >= 2 {
                        // ── 1. Lookahead ──────────────────────────────────────
                        // Only fires when EVERY arm provides all LOOKAHEAD_K bases,
                        // so score_arm_prefix compares the same number of bases on
                        // each arm.  Gating on the longest arm (as opposed to the
                        // shortest) let a short arm's score — capped at its own
                        // length — lose to a long arm's score purely because it
                        // was compared over fewer bases, regardless of which arm
                        // actually matched the query better.  A short arm (e.g. a
                        // length-1 "no insertion" shortcut against a many-node
                        // insertion arm) always lost this way, even when it was
                        // the correct one.  Short arms (< K) fall through to
                        // slide-and-lock's exhaustion rule or windowed DP, both of
                        // which score fairly.
                        let min_scored = all_arms
                            .iter()
                            .map(|arm| {
                                arm.len()
                                    .min(LOOKAHEAD_K)
                                    .min(query.len().saturating_sub(j_entry))
                            })
                            .min()
                            .unwrap_or(0);

                        let winner = if min_scored >= LOOKAHEAD_K {
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
                                    scores.iter().position(|&s| s == best)
                                } else {
                                    None
                                }
                            } else {
                                None
                            }
                        } else {
                            None
                        };

                        // ── 2. Slide-and-lock ─────────────────────────────────
                        let winner =
                            winner.or_else(|| slide_lock(&all_arms, nodes, query, j_entry));

                        // ── Mark losers + lock winning arm windows ────────────
                        if let Some(w) = winner {
                            for (arm_idx, arm) in all_arms.iter().enumerate() {
                                if arm_idx != w {
                                    for &losing_node in arm {
                                        lookahead_skip[losing_node] = true;
                                    }
                                }
                            }
                            // Give winning arm nodes tight per-depth j-windows
                            // so cells at the correct query column are always
                            // filled even when arm depth > spine_margin.
                            let arm = &all_arms[w];
                            let arm_len = arm.len();
                            for (d, &win_node) in arm.iter().enumerate() {
                                let j_c = (j_entry + d + 1).min(l) as u32;
                                scratch.lock_node_j.push((win_node, j_c));
                            }
                            scratch.lock_node_j.sort_unstable_by_key(|&(idx, _)| idx);
                            // Exit node: centred at the arm's expected terminal j.
                            let exit_node = topo[exit_t];
                            let exit_j = (j_entry + arm_len).min(l) as u32;
                            // Insert sorted (exit_node may already be present from
                            // a previous bubble in the same align() call).
                            match scratch
                                .lock_exit_j
                                .binary_search_by_key(&exit_node, |&(idx, _)| idx)
                            {
                                Ok(pos) => {
                                    scratch.lock_exit_j[pos].1 =
                                        exit_j.max(scratch.lock_exit_j[pos].1)
                                }
                                Err(pos) => scratch.lock_exit_j.insert(pos, (exit_node, exit_j)),
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
            let jlo = j_lo_arr[t];
            let jhi = j_hi_arr[t];
            let sm = gs(&m, t, l, jlo, jhi, row_width);
            let si = gs(&ins, t, l, jlo, jhi, row_width);
            let sd = gsd(&del, &del0, t, l, jlo, jhi, row_width);
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
        None => {
            // No node's banded row-window reached query column `l` at all --
            // the DP could not complete within the configured band. This
            // used to fall back to `(0, State::M)`, a fabricated starting
            // point nowhere near a real cell, which made the traceback
            // below immediately hit an UNSET cell and silently return an
            // empty `Vec<AlignOp>` -- no error, no warning, just a no-op
            // alignment that `add_to_graph` then applies as if the read
            // contributed nothing at all. Confirmed on real data
            // (`tests/adaptive_band_collapse.rs`): every non-seed read in a
            // population wider than the effective band hits this
            // identically, leaving the graph as just the seed's own linear
            // chain, which `consensus()`'s boundary trim then collapses to
            // a single base.
            //
            // `required` is estimated from `best_j_per_t`, which every row
            // already updates with the furthest query column its own
            // window reached (used elsewhere for the diagonal-skip
            // machinery, reused here rather than tracked separately): the
            // row that got closest to `l` fell short by
            // `l - max(best_j_per_t)`, so widening the margin by that same
            // amount is the natural single-shot estimate of what would
            // have worked. This is a heuristic, not a proof of
            // sufficiency -- see `align_with_retry`'s 3-pass caller, which
            // does not trust this estimate alone and escalates to a
            // mathematically-guaranteed-sufficient fully unbanded retry
            // if it turns out to still be too narrow.
            let max_best_j = best_j_per_t.iter().copied().max().unwrap_or(0);
            let shortfall = l.saturating_sub(max_best_j);
            let required = spine_margin.saturating_add(shortfall);
            return Err(PoaError::BandTooNarrow {
                configured: spine_margin,
                required,
            });
        }
    };

    // ── Traceback ─────────────────────────────────────────────────────────────
    let mut ops: Vec<AlignOp> = Vec::with_capacity(l + n / 4);
    let mut t = best_t;
    let mut j = l;
    let mut cur_state = best_state;

    // Defensive bound: a traceback over a well-formed DAG can never produce
    // more than `l` Insert/Match steps (each consumes one query base) plus
    // `n` Match/Delete steps (each consumes one graph node), so `l + n` is a
    // hard ceiling for any *correct* alignment. This crate has a documented
    // history of node-graph self-loops corrupting exactly this kind of
    // predecessor-pointer walk (see `design/graph_data_model_rework.md`
    // Phase 3's account of a literal node self-loop hanging
    // `heaviest_path`'s own traceback) -- if that class of bug is ever
    // reached from this traceback too, looping forever while `ops` grows
    // without bound is an OOM/crash risk, not just a slow answer. Treating
    // an over-long traceback as `BandTooNarrow` (rather than a silent
    // infinite loop) is conservative: it costs nothing in the overwhelming
    // majority of calls that never approach this bound, and turns a
    // process-ending crash into a catchable, retriable error for the rare
    // case that does.
    let max_ops = l.saturating_add(n).saturating_add(16);

    loop {
        if ops.len() > max_ops {
            return Err(PoaError::BandTooNarrow {
                configured: spine_margin,
                required: spine_margin.saturating_mul(2).max(l),
            });
        }
        let cell = {
            let jlo = j_lo_arr[t];
            let jhi = j_hi_arr[t];
            match cur_state {
                State::M => {
                    if j < jlo || j > jhi {
                        Cell::unset()
                    } else {
                        m[t * row_width + (j - jlo)]
                    }
                }
                State::I => {
                    if j < jlo || j > jhi {
                        Cell::unset()
                    } else {
                        ins[t * row_width + (j - jlo)]
                    }
                }
                State::D => {
                    if j == 0 {
                        del0[t]
                    } else if j < jlo || j > jhi {
                        Cell::unset()
                    } else {
                        del[t * row_width + (j - jlo)]
                    }
                }
            }
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
                t = cell.pred_t as usize;
                j -= 1;
                cur_state = best_prev_state_banded(
                    &m, &ins, &del, &del0, t, j, &j_lo_arr, &j_hi_arr, row_width,
                );
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
                let sm = gs(&m, t, j, j_lo_arr[t], j_hi_arr[t], row_width);
                let si = gs(&ins, t, j, j_lo_arr[t], j_hi_arr[t], row_width);
                cur_state = if sm >= si { State::M } else { State::I };
            }
            State::D => {
                ops.push(AlignOp::Delete(topo[t]));
                if cell.pred_t == VIRTUAL {
                    break;
                }
                t = cell.pred_t as usize;
                cur_state = best_prev_state_banded(
                    &m, &ins, &del, &del0, t, j, &j_lo_arr, &j_hi_arr, row_width,
                );
            }
        }

        if j == 0 && cur_state != State::D {
            break;
        }
    }

    ops.reverse();
    Ok(ops)
}

/// Wraps [`align`] with a bounded, guaranteed-terminating retry for
/// [`PoaError::BandTooNarrow`].
///
/// `TODO.md` has long described a "smart retry (3-pass)" here -- pass 1 at
/// the configured width, pass 2 at the error's own `required` estimate,
/// pass 3 fully unbanded -- as already implemented. It was not: `align`
/// never returned `Err` at all before the fix that added the
/// `BandTooNarrow` path above (confirmed by exhaustive search), so no
/// caller ever exercised a retry, and none existed. Several existing
/// tests (`band_too_narrow_fallback_to_unbanded`,
/// `tracking_band_survives_phase_shift`) already asserted the *outcome*
/// this retry produces, and passed anyway, for the wrong reason -- not
/// because a retry recovered a wide-enough band, but because the missing
/// `BandTooNarrow` path meant `align` silently returned a degenerate
/// (possibly empty) alignment that the assertions happened not to notice
/// was wrong. This is that design, actually built:
///
/// - **Pass 1**: the caller's own config, unmodified. This is the
///   overwhelmingly common path and costs nothing extra when it succeeds.
/// - **Pass 2**: on `BandTooNarrow`, retry once with `band_width` set to
///   the error's own `required` estimate and `adaptive_band` forced off
///   (so `required` is the literal effective margin used, not further
///   changed by the adaptive formula), and with `anchors` cleared -- an
///   anchor chain computed for the original, too-narrow window is not
///   necessarily valid for a wider one, and `align_read_ops`/
///   `align_read_ops_unbanded` already always call with no anchors for
///   exactly this reason (a "just get a correct alignment" call, not a
///   performance-sensitive incremental one).
/// - **Pass 3**: if pass 2 *also* returns `BandTooNarrow` (the pass-1
///   heuristic `required` estimate under-shot -- `real_in_edge_count`-style
///   estimates elsewhere in this file are deliberately conservative, but
///   this one is a single-shot approximation, not a proven lower bound),
///   retry with a fully unbanded config (`band_width = 0, adaptive_band =
///   false`, i.e. `spine_margin = l`) -- mathematically guaranteed
///   sufficient (full NW over the DAG), not another estimate. If this
///   *still* errors, that is a genuine, unexpected failure and is
///   propagated rather than looped on again.
///
/// Returns `(ops, retried)`: `retried` is `true` whenever pass 1 (the
/// caller's own config) did not succeed on its own, so the caller can tell
/// a clean single-band alignment from one that only succeeded after
/// widening -- see `PoaGraph::used_band_retry` for why that distinction
/// matters beyond just "did it eventually succeed."
#[allow(clippy::too_many_arguments)]
fn align_with_retry(
    nodes: &[Node],
    topo: &[usize],
    rank_of: &[usize],
    spine: &[(usize, u8, i32)],
    query: &[u8],
    cfg: &PoaConfig,
    scratch: &mut AlignScratch,
    anchors: &[(usize, usize)],
) -> Result<(Vec<AlignOp>, bool), PoaError> {
    match align(nodes, topo, rank_of, spine, query, cfg, scratch, anchors) {
        Ok(ops) => Ok((ops, false)),
        Err(PoaError::BandTooNarrow { required, .. }) => {
            let mut cfg2 = cfg.clone();
            cfg2.band_width = required;
            cfg2.adaptive_band = false;
            match align(nodes, topo, rank_of, spine, query, &cfg2, scratch, &[]) {
                Ok(ops) => Ok((ops, true)),
                Err(PoaError::BandTooNarrow { .. }) => {
                    let mut cfg3 = cfg.clone();
                    cfg3.band_width = 0;
                    cfg3.adaptive_band = false;
                    align(nodes, topo, rank_of, spine, query, &cfg3, scratch, &[])
                        .map(|ops| (ops, true))
                }
                Err(other) => Err(other),
            }
        }
        Err(other) => Err(other),
    }
}

#[allow(clippy::too_many_arguments)]
#[inline]
fn best_prev_state_banded(
    m: &[Cell],
    ins: &[Cell],
    del: &[Cell],
    del0: &[Cell],
    t: usize,
    j: usize,
    j_lo_arr: &[usize],
    j_hi_arr: &[usize],
    row_width: usize,
) -> State {
    let jlo = j_lo_arr[t];
    let jhi = j_hi_arr[t];
    let sm = gs(m, t, j, jlo, jhi, row_width);
    let si = gs(ins, t, j, jlo, jhi, row_width);
    let sd = gsd(del, del0, t, j, jlo, jhi, row_width);
    if sm != UNSET && sm >= si && sm >= sd {
        State::M
    } else if si != UNSET && si >= sd {
        State::I
    } else {
        State::D
    }
}

// ─── Graph update ─────────────────────────────────────────────────────────────

/// Verifies (read-only) that the existing chain recorded at `start` still
/// spells exactly `edit`: every one of the `edit.len() - 1` remaining hops
/// must follow a single out-edge (no internal branching since the chain was
/// recorded) to a node whose base matches the corresponding byte of `edit`.
/// Returns the full chain of node indices on success, or `None` if the chain
/// has since branched or drifted -- the caller falls back to creating a
/// fresh arm exactly as before this phase existed.
fn verify_reuse_chain(nodes: &[Node], edit: &[u8], start: usize) -> Option<Vec<usize>> {
    debug_assert!(!edit.is_empty());
    let mut chain = Vec::with_capacity(edit.len());
    let mut cur = start;
    if nodes[cur].base != edit[0] {
        return None;
    }
    chain.push(cur);
    for &b in &edit[1..] {
        match nodes[cur].out_edges.as_slice() {
            [(next, _)] if nodes[*next].base == b => {
                cur = *next;
                chain.push(cur);
            }
            _ => return None,
        }
    }
    Some(chain)
}

/// Returns `true` if any node in `chain` is referenced again anywhere in
/// `rest` (as a `Match` or `Delete` target).
///
/// This is the safety check that closes the gap the naive first attempt at
/// this feature missed and this phase's own first implementation attempt
/// also missed: `align()` computes the *entire* traceback for a read in one
/// pass, before `add_to_graph` runs at all, and its `Match`/`Delete` ops
/// reference specific existing node indices that the DP decided this read
/// visits -- independently of whatever `add_to_graph` does. If content-address
/// reuse redirects `prev` to an existing node, and that *same* node is also
/// named later in this read's own traceback (which was computed without any
/// knowledge that reuse would happen), committing the reuse creates a second,
/// unrelated edge into that node from `prev` -- up to and including a literal
/// self-loop when the very next op names the node reuse just landed on.
/// Confirmed concretely on real data (DAB1 SCA37, seed=25): node 49 ended up
/// with `in_edges=[17, 49]` and `out_edges=[18, 99, 49]`, a genuine self-loop,
/// which corrupts `topological_order`'s Kahn's-algorithm bookkeeping (in-degree
/// can never reach zero) and hangs `heaviest_path`'s traceback (an unbounded
/// walk over a predecessor-pointer array that develops its own cycle once
/// `rank_of` assigns the same rank to two different nodes). Skipping reuse
/// whenever the candidate chain reappears later in this same read's ops is a
/// cheap (bounded by `ops.len()`), always-correct guard against exactly this.
fn reuse_would_collide(chain: &[usize], rest: &[AlignOp]) -> bool {
    rest.iter().any(|op| match *op {
        AlignOp::Match(idx) | AlignOp::Delete(idx) => chain.contains(&idx),
        AlignOp::Insert(_) => false,
    })
}

/// Returns `true` if committing this reuse would create a back edge relative
/// to the *rest of this same read's own traceback*, not just relative to the
/// fork itself (that narrower check lives in `try_reuse_arm`, comparing only
/// `fork` against `chain[0]`).
///
/// `align()`'s own traceback is provably rank-monotonic in the *original*
/// node identities it names: reading its `Match`/`Delete` ops in forward
/// order, the referenced existing-node rank strictly increases throughout
/// (each step in the backward-built traceback moves to a strictly smaller
/// DP row, so after `ops.reverse()` the forward sequence strictly
/// increases). Confirmed empirically by dumping a real failing read's ops
/// end to end: no adjacent pair of `Match`/`Delete` references ever regresses
/// in rank.
///
/// Reuse breaks this guarantee from the *outside*: it substitutes a
/// *different*, already-existing node (`chain`'s last element, at whatever
/// rank that node happens to have) in place of whatever fresh identity
/// `align()`'s own traceback implicitly expected at this position. The very
/// next `Match`/`Delete` op in `rest` (an ordinary, unprotected
/// `increment_or_add_edge` call in `add_to_graph`, ordinary edges never
/// having needed a rank check before reuse existed) then wires an edge from
/// this substitute straight to whatever existing node `align()` visits next
/// -- which is only ever safe if the substitute's rank is still *less than*
/// that next node's rank. `chain`'s internal edges are all real, pre-existing
/// graph edges (verified by `verify_reuse_chain`'s single-out-edge walk), so
/// rank strictly increases along `chain` too; checking only the *last*
/// element (the one `prev` becomes) against the first subsequent existing-
/// node reference is therefore sufficient -- if that holds, every earlier
/// `chain` element (smaller rank) holds too.
///
/// Confirmed concretely via a `cag50_hifi`-class read set once the
/// `BandTooNarrow` retry started actually widening bands for some reads: a
/// substitution reuse landed on an existing node whose rank happened to sit
/// *after* the very next node `align()`'s own traceback visited, splicing in
/// a real back edge and corrupting the graph into a cycle -- silently
/// hanging `heaviest_path`'s traceback in an unbounded loop (the same
/// failure class `reuse_would_collide` already guards, reached by a
/// different route: rank inversion rather than direct index reappearance).
fn reuse_would_create_back_edge(rank_of: &[usize], chain: &[usize], rest: &[AlignOp]) -> bool {
    let Some(&last_rank) = chain.last().and_then(|&idx| rank_of.get(idx)) else {
        return false;
    };
    for op in rest {
        if let AlignOp::Match(idx) | AlignOp::Delete(idx) = *op {
            return match rank_of.get(idx) {
                Some(&r) => last_rank >= r,
                // `idx` isn't in this read's entry-snapshot rank_of at all,
                // meaning it's a node this same read created after the
                // snapshot was taken -- Match/Delete ops from align()'s own
                // traceback never reference such nodes (they only name
                // nodes that existed before this read started), so this
                // shouldn't happen; bail out defensively rather than reject.
                None => false,
            };
        }
    }
    false
}

/// Commits a verified reuse chain: increments edge weight/`edge_reads` and
/// bumps `coverage` on every reused node -- the same effect a genuine `Match`
/// against each of these nodes would have had. Only ever called after both
/// `verify_reuse_chain` and `reuse_would_collide` have cleared the chain.
fn commit_reuse_chain(
    nodes: &mut [Node],
    edge_reads: &mut HashMap<(usize, usize), Vec<u32>>,
    fork: usize,
    chain: &[usize],
    read_idx: u32,
) -> usize {
    let mut prev_node = fork;
    for &node in chain {
        // Reuse only ever traverses an edge that `verify_reuse_chain` already
        // confirmed exists (that's what makes it a "reuse" and not a fresh
        // node creation) -- this call should always increment, never create.
        // A brand-new edge appearing here would mean a reused arm's fork
        // propagation was never run when the arm was first built, which
        // would be a real bug elsewhere, not something to paper over here.
        let created = increment_or_add_edge(nodes, prev_node, node, false);
        debug_assert!(
            !created,
            "reuse chain traversed an edge that did not already exist"
        );
        edge_reads
            .entry((prev_node, node))
            .or_default()
            .push(read_idx);
        nodes[node].coverage += 1;
        prev_node = node;
    }
    *chain.last().unwrap()
}

/// Looks up, verifies, and safety-checks a candidate reuse in one call;
/// returns the reused chain's last node index on success. See
/// `verify_reuse_chain` and `reuse_would_collide` for what "success" means.
///
/// `rank_of` is the topological rank snapshot taken at the *start* of the
/// current `add_read` call (before this read's own ops are applied). Reuse
/// splices a brand-new edge `fork -> chain[0]` between two nodes that both
/// already existed in the graph -- unlike an ordinary Match/Insert/Delete op,
/// which only ever connects nodes along `align()`'s own traceback and is
/// therefore guaranteed rank-monotonic by construction (each step moves to a
/// strictly smaller DP row), this edge is *not* derived from that traceback
/// at all. It comes from `fork_arm_index`, a lookup keyed only by content
/// (the fork node and the edit's bytes), with no reference to where either
/// node sits in the graph's real topological order relative to the other.
///
/// `reuse_would_collide` guards one specific way this can go wrong (the
/// reused chain reappearing later in this same read's own remaining ops),
/// but not this one: if `fork` is not actually an ancestor of `chain[0]` in
/// the graph's existing structure (i.e. `chain[0]` already sits at or before
/// `fork` in topological order), splicing in `fork -> chain[0]` creates a
/// genuine multi-node cycle, not a same-read self-loop. `reuse_would_collide`
/// cannot see this because it only ever looks within the current read's own
/// ops; it has no way to know the reused chain's actual position in the rest
/// of the graph.
///
/// Confirmed concretely via a `cag50_hifi`-class read set once the
/// `BandTooNarrow` retry (see `align_with_retry`) started actually firing:
/// a retried alignment pass (wider band, anchors cleared) can settle on a
/// different fork/traceback shape than the pass that originally populated
/// `fork_arm_index` for this position, so by the time a later read's ops
/// reach this same fork, `chain[0]` may already be topologically upstream of
/// `fork` rather than downstream. Committing the reuse in that case corrupts
/// the graph into a real cycle: `topological_order`'s Kahn's-algorithm queue
/// never reaches in-degree zero for the cyclic nodes, silently excluding them
/// from `topo` (and leaving their `rank_of` entries at the default `0`,
/// colliding with whatever node is legitimately first), which in turn hangs
/// `heaviest_path`'s own predecessor-pointer traceback in an unbounded loop
/// -- the same failure class already documented on `reuse_would_collide`,
/// just reached by a different path. Rejecting the reuse whenever it would
/// not be rank-increasing is cheap (an O(1) rank comparison) and always
/// safe: falling through to fresh-node creation, exactly like the existing
/// `None`-returning cases, can never introduce a cycle since a brand-new node
/// has no edges to anywhere except `fork` and, later, its own single
/// successor.
#[allow(clippy::too_many_arguments)]
fn try_reuse_arm(
    nodes: &mut [Node],
    edge_reads: &mut HashMap<(usize, usize), Vec<u32>>,
    fork_arm_index: &HashMap<usize, HashMap<Vec<u8>, usize>>,
    rank_of: &[usize],
    fork: usize,
    edit: &[u8],
    rest_ops: &[AlignOp],
    read_idx: u32,
) -> Option<usize> {
    let start = *fork_arm_index.get(&fork)?.get(edit)?;
    if fork >= rank_of.len() || start >= rank_of.len() || rank_of[fork] >= rank_of[start] {
        return None;
    }
    let chain = verify_reuse_chain(nodes, edit, start)?;
    if reuse_would_collide(&chain, rest_ops) {
        return None;
    }
    if reuse_would_create_back_edge(rank_of, &chain, rest_ops) {
        return None;
    }
    Some(commit_reuse_chain(
        nodes, edge_reads, fork, &chain, read_idx,
    ))
}

/// Record (or increment) a bypass edge `from -> to`. See
/// `PoaGraph.bypass_edges`. Weight is a running count of reads that took this
/// bypass; a self-loop (`from == to`) is refused defensively for symmetry
/// with `increment_or_add_edge`, though the caller's own guards should make
/// it unreachable.
fn record_bypass_edge(
    bypass_edges: &mut HashMap<usize, Vec<(usize, i32)>>,
    from: usize,
    to: usize,
) {
    if from == to {
        debug_assert_ne!(from, to, "record_bypass_edge: refusing a self-loop bypass");
        return;
    }
    let entry = bypass_edges.entry(from).or_default();
    for (t, w) in entry.iter_mut() {
        if *t == to {
            *w += 1;
            return;
        }
    }
    entry.push((to, 1));
}

// `edge_delete_reads` is intentionally NOT a parameter: under pure bypass a
// deleting read touches the skipped node only via `delete_count`, never via an
// edge, so there is nothing to record in that (write-only, unread) map. The
// `PoaGraph.edge_delete_reads` field is left in place (initialised empty,
// never populated) as a removable-later cleanup, out of scope here.
#[allow(clippy::too_many_arguments)]
fn add_to_graph(
    nodes: &mut Vec<Node>,
    edge_reads: &mut HashMap<(usize, usize), Vec<u32>>,
    bypass_edges: &mut HashMap<usize, Vec<(usize, i32)>>,
    fork_arm_index: &mut HashMap<usize, HashMap<Vec<u8>, usize>>,
    rank_of: &[usize],
    query: &[u8],
    ops: &[AlignOp],
    read_idx: u32,
) {
    let mut prev: Option<usize> = None;
    let mut q_idx: usize = 0;
    let mut i = 0usize;

    // Bypass-edge tracking (Phase 1 of design/bypass_edge_delete_rework.md).
    // `bypass_pending` is `Some(entry_pred)` while inside a run of consecutive
    // `Delete` ops, where `entry_pred` is `prev` as it stood *just before* the
    // run's first Delete (itself an `Option`: `None` when the read begins with
    // a Delete run, i.e. there is no predecessor to bypass *from*). It is set
    // on the first Delete of a run, left untouched by subsequent Deletes in
    // the same run, and consumed by the next `Match`/`Insert` op -- which, by
    // the time it finishes, has advanced `prev` to the node the read resumed
    // at (the bypass edge's `to`). Capturing the resume node this way, rather
    // than peeking at `ops[i]`, is what makes it correct in the mismatch and
    // insert cases too: those create/reuse a *new* node whose index isn't
    // known until the op is processed, but `prev` names it correctly
    // afterward regardless. A run left pending at end-of-ops (a terminal
    // Delete run, e.g. a semi-global read that runs out before the graph does)
    // is simply dropped -- there is no resume node, so no bypass edge to
    // nowhere is created.
    let mut bypass_pending: Option<Option<usize>> = None;

    // Nodes that gained a new out-edge while processing *this* read.
    // `propagate_fork_if_new` is deliberately NOT called inline as each edge
    // is added -- confirmed empirically (see `set_new_node_own_fork`'s doc
    // comment) that doing so races against this same read's own later ops:
    // a node can look like an ordinary single-predecessor node at the moment
    // a fork upstream is created, then gain a second in-edge (reconvergence)
    // a few ops later in this exact read's own traceback. Collecting
    // candidates here and propagating once, after the whole read's ops are
    // processed, means `in_edges.len()` already reflects this read's final
    // effect on the graph by the time any propagation walk inspects it.
    let mut newly_forked: Vec<usize> = Vec::new();

    while i < ops.len() {
        // The node this iteration's op reconnects the read to (its first, in
        // read order -- not `prev`'s final value, which for a multi-base
        // Insert is the run's *last* node). Set (to `Some`) by the
        // `Match`/`Insert` arms to signal "a non-Delete op ran, close any
        // pending bypass run"; left `None` by `Delete` so a run stays open
        // across its whole span.
        let mut resume_node: Option<usize> = None;
        // Whether this iteration's op was a clean Match onto a node ALREADY in
        // the graph (base agrees) -- the only resume shape that records a
        // bypass edge under pure bypass. A mismatch/insert (new or reused
        // structure) instead keeps its ordinary real out-edge and records no
        // bypass (see the resolution block after the match).
        let mut resume_is_bypass = false;
        match ops[i] {
            AlignOp::Match(node_idx) => {
                let q_base = query[q_idx];
                q_idx += 1;
                let cur = if nodes[node_idx].base == q_base {
                    nodes[node_idx].coverage += 1;
                    if bypass_pending.is_some() {
                        // Pure bypass (revised Phase 1): this clean Match is
                        // resuming a Delete run onto a node already in the
                        // graph. Record NO matched out-edge (from the entry
                        // predecessor or anywhere) and NO `edge_reads` entry --
                        // that edge is exactly the "laundering" this rework
                        // removes (it let the skipped node's downstream matched
                        // weight re-inflate its arm). `coverage` is still bumped
                        // (above); the reconnection is represented solely by the
                        // bypass edge recorded at the resolution block below.
                        resume_is_bypass = true;
                    } else if let Some(p) = prev {
                        if increment_or_add_edge(nodes, p, node_idx, false) {
                            newly_forked.push(p);
                        }
                        edge_reads.entry((p, node_idx)).or_default().push(read_idx);
                    }
                    node_idx
                } else {
                    // Single-base substitution edit. Try exact-duplicate reuse
                    // at this fork before creating a new node.
                    let edit = [q_base];
                    let reused = prev.and_then(|p| {
                        try_reuse_arm(
                            nodes,
                            edge_reads,
                            fork_arm_index,
                            rank_of,
                            p,
                            &edit,
                            &ops[i + 1..],
                            read_idx,
                        )
                    });
                    if let Some(reused_idx) = reused {
                        reused_idx
                    } else {
                        let new_idx = push_node(nodes, q_base);
                        nodes[new_idx].coverage = 1;
                        if let Some(p) = prev {
                            nodes[p].out_edges.push((
                                new_idx,
                                EdgeWeight {
                                    matched: 1,
                                    deleted: 0,
                                },
                            ));
                            nodes[new_idx].in_edges.push(p);
                            edge_reads.entry((p, new_idx)).or_default().push(read_idx);
                            fork_arm_index
                                .entry(p)
                                .or_default()
                                .insert(edit.to_vec(), new_idx);
                            set_new_node_own_fork(nodes, p, new_idx);
                            newly_forked.push(p);
                        }
                        new_idx
                    }
                };
                // A Match reconnects at exactly one node (`cur`), whether a
                // clean match, a reused single-base substitution (chain
                // length 1, so first == last), or a freshly created
                // substitute -- so `cur` is unambiguously the resume node.
                resume_node = Some(cur);
                prev = Some(cur);
                i += 1;
            }
            AlignOp::Insert(first_base) => {
                // Collect the full run of consecutive Insert ops: the
                // "characterized edit" is only fully known once the run
                // ends, so the content-address lookup/insert has to happen
                // for the whole run at once, not per base.
                let run_start = i;
                let mut edit = vec![first_base];
                let mut j = i + 1;
                while j < ops.len() {
                    if let AlignOp::Insert(b) = ops[j] {
                        edit.push(b);
                        j += 1;
                    } else {
                        break;
                    }
                }

                let reused = prev.and_then(|p| {
                    try_reuse_arm(
                        nodes,
                        edge_reads,
                        fork_arm_index,
                        rank_of,
                        p,
                        &edit,
                        &ops[j..],
                        read_idx,
                    )
                });
                if let Some(reused_idx) = reused {
                    // An Insert creates/rejoins DIVERGENT structure, not the
                    // existing main path, so under pure bypass it keeps its
                    // ordinary real edges (created inside `try_reuse_arm`) and
                    // records NO bypass edge -- `resume_is_bypass` stays false.
                    // `resume_node` just needs to be `Some` to signal that a
                    // non-Delete op ran and any pending bypass run should close.
                    resume_node = Some(reused_idx);
                    prev = Some(reused_idx);
                } else {
                    let fork = prev;
                    let mut chain_start = None;
                    for &b in &edit {
                        let new_idx = push_node(nodes, b);
                        nodes[new_idx].coverage = 1;
                        if let Some(p) = prev {
                            nodes[p].out_edges.push((
                                new_idx,
                                EdgeWeight {
                                    matched: 1,
                                    deleted: 0,
                                },
                            ));
                            nodes[new_idx].in_edges.push(p);
                            edge_reads.entry((p, new_idx)).or_default().push(read_idx);
                            set_new_node_own_fork(nodes, p, new_idx);
                            newly_forked.push(p);
                        }
                        chain_start.get_or_insert(new_idx);
                        prev = Some(new_idx);
                    }
                    if let (Some(f), Some(start)) = (fork, chain_start) {
                        fork_arm_index
                            .entry(f)
                            .or_default()
                            .insert(edit.clone(), start);
                    }
                    // New structure, not the existing main path: keep the real
                    // edges just created, record no bypass (`resume_is_bypass`
                    // stays false); `resume_node` only signals run closure.
                    resume_node = chain_start;
                }
                q_idx += edit.len();
                i = run_start + edit.len();
                let _ = j; // j == i; kept named for clarity above
            }
            AlignOp::Delete(node_idx) => {
                // Open a bypass run on the first Delete: capture `prev` (the
                // entry predecessor). `Some(None)` marks a run that begins with
                // no predecessor (read starts on a Delete under semi-global),
                // which resolves to no bypass edge.
                if bypass_pending.is_none() {
                    bypass_pending = Some(prev);
                }
                // Pure bypass (revised Phase 1): increment the skipped node's
                // own `delete_count` (this is all `mean_column_entropy` and
                // `majority_frequency` need -- they read per-node
                // coverage/delete_count only) and do NOTHING else. In
                // particular do NOT:
                //   - create/increment a `p -> node_idx` edge
                //     (`increment_or_add_edge(.., true)`), nor push
                //     `edge_delete_reads` (that map is write-only -- no
                //     production consumer reads it -- so under pure bypass it
                //     simply stops being populated; the field is left in place,
                //     removable in a later cleanup);
                //   - advance `prev` through the skipped node.
                // Leaving `prev` at the entry predecessor is what makes the
                // read's next real op reconnect via a bypass edge (existing
                // resume) or a fresh minority out-edge (new resume) instead of
                // laundering the skipped node's downstream matched edge.
                nodes[node_idx].delete_count += 1;
                i += 1;
            }
        }

        // Resolve a pending bypass run at the first non-Delete op. A `Delete`
        // leaves `resume_node` `None`, so the run stays open across its whole
        // span and is only closed here by the Match/Insert that ends it. A
        // bypass edge is recorded only when that resume was a clean Match onto
        // an EXISTING node (`resume_is_bypass`) AND the run had an entry
        // predecessor to bypass *from* (a leading Delete run has none). A
        // mismatch/insert resume closes the run but records no bypass -- its
        // real out-edge already represents the reconnection.
        if bypass_pending.is_some() && resume_node.is_some() {
            if resume_is_bypass {
                if let (Some(Some(from)), Some(to)) = (bypass_pending, resume_node) {
                    // `rank_of` is the pre-read topological snapshot. A bypass
                    // run whose entry predecessor `from` is a node this SAME
                    // read created earlier (e.g. an Insert that ran before the
                    // Delete run) has `from >= rank_of.len()` and cannot be
                    // rank-checked against the snapshot. Such a bypass is
                    // forward by construction -- ops are applied in alignment
                    // order, so an Insert-created `from` precedes the
                    // later-matched existing `to` -- so only assert the
                    // topological order when both endpoints are in the
                    // snapshot. Mirrors `try_reuse_arm`'s own `>= rank_of.len()`
                    // guard. Debug-only: the recorded bypass edge is identical
                    // either way (this was a spurious `debug_assert` panic in
                    // debug/CI builds; release builds, where `debug_assert` is
                    // compiled out, always recorded the correct edge).
                    debug_assert!(
                        from >= rank_of.len() || to >= rank_of.len() || rank_of[from] < rank_of[to],
                        "bypass edge {from}->{to} between in-snapshot nodes must respect \
                         topological order -- it is redundant with the real \
                         entry-pred->...->resume path through the skipped nodes"
                    );
                    // Anti-laundering guard: a bypass resume must NOT have
                    // created a matched edge into the resume node for this
                    // read. If it had, this read would appear as the last
                    // pushed reader on that edge -- the exact laundering the
                    // pure-bypass Delete arm removes.
                    debug_assert!(
                        edge_reads
                            .get(&(from, to))
                            .is_none_or(|v| v.last() != Some(&read_idx)),
                        "laundering guard: deleting read {read_idx} created a matched \
                         edge {from}->{to} into its bypass resume node"
                    );
                    record_bypass_edge(bypass_edges, from, to);
                }
            }
            bypass_pending = None;
        }
    }

    // A bypass run still open here is a terminal Delete run (the read ran out
    // of query before the graph ended, e.g. a semi-global right-flank gap):
    // no resume node ever followed, so there is no bypass edge to record.
    // Dropping `bypass_pending` unresolved is the correct handling -- a bypass
    // edge to a nonexistent "next" node would be malformed.
    let _ = bypass_pending;

    // Now that this read's entire traceback has been applied and the graph
    // reflects its final shape (including any reconvergence this same read
    // itself created), it's safe to propagate fork-context to descendants.
    // Dedup first: a fork can appear multiple times (e.g. an Insert run's
    // per-base loop calls this for the same `p` only once, but a read could
    // still revisit the same fork node as `p` more than once via separate
    // ops), and `propagate_fork_if_new` is cheap but there's no reason to
    // redo the same walk twice.
    newly_forked.sort_unstable();
    newly_forked.dedup();
    for p in newly_forked {
        propagate_fork_if_new(nodes, p);
    }
}

// ─── Heaviest path ────────────────────────────────────────────────────────────

/// Cumulative-weight DP over the DAG, picking the heaviest incoming edge at
/// each node to build the consensus spine.
///
/// Scores on `edge.matched` only -- `deleted` traversals contribute nothing,
/// not even a discount. A read that deletes through a node provides zero
/// evidence that the node's *base* is correct (it skipped confirming it
/// entirely); scoring on total (matched + deleted) traffic let a node reached
/// mostly by Delete out-compete a genuinely-Match-confirmed alternative arm
/// at the same fork, confirmed on a forced case (6 reads deleting through a
/// reference base vs. 4 reads genuinely matching a SNP base at the same
/// position -- see `heaviest_path_prefers_matched_over_delete_inflated_arm`
/// in src/tests.rs). This mirrors the interior filter's own existing
/// philosophy (`coverage > delete_count`): Delete already counts as zero
/// evidence for node inclusion there, so it would be inconsistent for this
/// function's arm-*selection* to treat it as partial evidence instead.
///
/// **Bypass edges (Phase 2 of `design/bypass_edge_delete_rework.md`) are the
/// mechanism that lets the DP route around a widely-skipped node.** A run of
/// reads that Deleted past a node reconnected via a `bypass_edges` entry
/// `from -> to` (recorded in Phase 1); here each such entry is relaxed as an
/// additional competing outgoing option of `from`, scored `(weight - 1)` on
/// the plain `i32` bypass weight directly -- NOT via `EdgeWeight.matched`, a
/// bypass weight does not live in an `EdgeWeight` at all. So a node most reads
/// skip is out-competed by the bypass around it, exactly as abPOA/SPOA's plain
/// heaviest-bundling does, and there is no laundering to counteract: under pure
/// bypass (Phase 1) the deleting reads never touch the skipped node's
/// downstream matched edge, so that edge cannot re-inflate the weak node's arm
/// at a reconvergence point. This is why the old `credibility_penalty` (which
/// discounted net-Delete-dominant nodes to counter exactly that laundering) is
/// **retired** here -- Audit item 1's redundancy proof holds once the laundered
/// edge is gone. (`heaviest_path_prefers_matched_over_delete_inflated_arm` is
/// the ground truth: it now lands on the SNP arm via the bypass edge + no
/// penalty, rather than via the penalty on a de-laundered graph.)
///
/// A bypass edge `(from, to)` always satisfies `rank_of[from] < rank_of[to]`
/// (it is redundant with the real `from -> ...skipped... -> to` path that
/// still exists in the graph, so it can introduce no cycle) -- asserted below.
/// The DP is forward-relaxation in topological order, so relaxing `from`'s
/// bypass edges into `cum[rank_of[to]]` when `from` is processed correctly
/// makes the bypass an incoming option for `to` (its `cum` is not finalized
/// until every strictly-earlier rank, `from` included, has been processed).
fn heaviest_path(
    nodes: &[Node],
    topo: &[usize],
    rank_of: &[usize],
    bypass_edges: &HashMap<usize, Vec<(usize, i32)>>,
) -> Vec<(usize, u8, i32)> {
    let n = topo.len();
    let mut cum: Vec<(i64, Option<usize>, i32)> = vec![(0, None, 0); n];

    for t in 0..n {
        let node_idx = topo[t];
        let node = &nodes[node_idx];
        let curr = cum[t].0;
        for &(succ_idx, ew) in &node.out_edges {
            let succ_t = rank_of[succ_idx];
            let candidate = curr + (ew.matched - 1) as i64;
            if candidate > cum[succ_t].0 {
                cum[succ_t] = (candidate, Some(t), ew.matched);
            }
        }
        // Bypass edges are relaxed AFTER out-edges, so an exact tie between a
        // bypass and a through-edge to the same target keeps the through-edge
        // (strict `>` below) -- i.e. prefer keeping a node over skipping it
        // when the evidence is exactly balanced, leaving such a node on the
        // spine for the interior filter to judge (pre-Phase-4 behaviour).
        if let Some(bypasses) = bypass_edges.get(&node_idx) {
            for &(succ_idx, weight) in bypasses {
                let succ_t = rank_of[succ_idx];
                debug_assert!(
                    t < succ_t,
                    "bypass edge {node_idx}->{succ_idx} (rank {t}->{succ_t}) must respect \
                     topological order; it is redundant with the real \
                     from->...->to path through the skipped nodes, so a violation means \
                     the graph is not a DAG"
                );
                let candidate = curr + (weight - 1) as i64;
                if candidate > cum[succ_t].0 {
                    cum[succ_t] = (candidate, Some(t), weight);
                }
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

    // Phase 2 territory (design/graph_data_model_rework.md): kept on total
    // (matched + deleted) weight here, unchanged from pre-split behaviour.
    // Splitting this by traversal type is deferred to the GraphStats/
    // find_bubbles rework, not part of this pass.
    let mut weights: Vec<f64> = nodes
        .iter()
        .flat_map(|nd| nd.out_edges.iter().map(|&(_, ew)| ew.total() as f64))
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
    let mut longest_bubble_span = 0usize;
    // Phase 2 territory: same total-weight threshold as before the split.
    for nd in nodes {
        let qualifying: Vec<(usize, i32)> = nd
            .out_edges
            .iter()
            .filter(|&&(_, ew)| ew.total() >= threshold)
            .map(|&(to, ew)| (to, ew.total()))
            .collect();
        if qualifying.len() >= 2 {
            bubble_count += 1;
            let mut weights: Vec<i32> = qualifying.iter().map(|&(_, w)| w).collect();
            weights.sort_unstable_by(|a, b| b.cmp(a));
            max_bubble_depth = max_bubble_depth.max(weights[1] as usize);
            for &(arm_start, _) in &qualifying {
                let span = if nodes[arm_start].in_edges.len() > 1 {
                    0 // direct edge to exit: 0-length arm
                } else {
                    materialize_arm_len(nodes, arm_start, ARM_MAX_DEPTH)
                };
                longest_bubble_span = longest_bubble_span.max(span);
            }
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
        longest_bubble_span,
        median_input_read_len: 0, // caller fills this in from self.reads
    }
}

fn median_read_len(reads: &[Vec<u8>]) -> usize {
    if reads.is_empty() {
        return 0;
    }
    let mut lens: Vec<usize> = reads.iter().map(|r| r.len()).collect();
    lens.sort_unstable();
    lens[lens.len() / 2]
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
    // Thresholded on matched-only weight: a delete-heavy pseudo-arm (reads
    // that skip past this position, not reads that confirm a different
    // base) should not qualify as a competing allele candidate.
    topo.iter()
        .copied()
        .filter_map(|node_idx| {
            let arms: Vec<usize> = nodes[node_idx]
                .out_edges
                .iter()
                .filter(|&&(_, ew)| ew.matched >= threshold)
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

/// Finds bubbles where at least one arm has structural size (span ≥ cfg.phasing_bubble_min_span).
/// These represent genuine length variants or SVs rather than SNPs / short indels.
/// Used for cross-bubble compatibility phasing in consensus_multi.
fn find_structural_bubbles(
    nodes: &[Node],
    topo: &[usize],
    n_reads: usize,
    cfg: &PoaConfig,
) -> Vec<(usize, Vec<usize>)> {
    let threshold = ((n_reads as f64 * cfg.min_allele_freq).ceil() as i32).max(1);
    let max_check = cfg.phasing_bubble_min_span.saturating_add(1);

    topo.iter()
        .copied()
        .filter_map(|entry_node| {
            // Matched-only, same reasoning as find_bubbles above.
            let arms: Vec<usize> = nodes[entry_node]
                .out_edges
                .iter()
                .filter(|&&(_, ew)| ew.matched >= threshold)
                .map(|&(to, _)| to)
                .collect();

            if arms.len() < 2 {
                return None;
            }

            // An arm whose start node already has a REAL (significant)
            // second in-edge is a direct edge to the bubble exit (arm span
            // 0); a start node with only NOISE extra in-edges is not
            // actually a reconvergence and its span is still measured
            // normally (see materialize_arm_len_tolerant's doc comment).
            let max_span = arms
                .iter()
                .map(|&start| {
                    if real_in_edge_count(nodes, start, threshold) > 1 {
                        0
                    } else {
                        materialize_arm_len_tolerant(nodes, start, max_check, threshold)
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

// ─── Bubble site extraction ───────────────────────────────────────────────────

/// Maximum arm sequence length captured in [`BubbleSite::arm_sequences`].
/// Arms longer than this return an empty `Vec`; callers use `arm_read_counts`
/// to assess support for large structural variants.
const ARM_SEQUENCE_CAP: usize = 256;

/// Collect the base sequence of one bubble arm from `start` up to (but not
/// including) the exit/reconvergence node.
///
/// Returns an empty `Vec` when:
/// - `start` is itself the exit (multiple in-edges → 0-length deletion arm).
/// - The arm exceeds `ARM_SEQUENCE_CAP` nodes.
/// - The arm has an internal branch (nested bubble).
fn arm_sequence(nodes: &[Node], start: usize) -> Vec<u8> {
    if nodes[start].in_edges.len() > 1 {
        return vec![];
    }
    let mut seq = Vec::new();
    let mut cur = start;
    loop {
        if seq.len() >= ARM_SEQUENCE_CAP {
            return vec![];
        }
        seq.push(nodes[cur].base);
        match nodes[cur].out_edges.as_slice() {
            [(next, _)] => {
                if nodes[*next].in_edges.len() > 1 {
                    break; // next is the exit node
                }
                cur = *next;
            }
            _ => break, // dead end or internal branch; return what we have
        }
    }
    seq
}

/// Build the [`BubbleSite`] list for a consensus path.
///
/// Scans the graph for nodes with two or more outgoing arms each meeting the
/// `min_allele_freq` threshold.  Only nodes that appear on `filtered` (the
/// trimmed consensus path) are included; off-path bubble entries are skipped.
fn collect_bubble_sites(
    nodes: &[Node],
    topo: &[usize],
    filtered: &[(usize, u8, i32)],
    edge_reads: &HashMap<(usize, usize), Vec<u32>>,
    min_allele_freq: f64,
    n_reads: usize,
    phasing_bubble_min_span: usize,
) -> Vec<BubbleSite> {
    let path_pos: HashMap<usize, usize> = filtered
        .iter()
        .enumerate()
        .map(|(i, &(node_idx, _, _))| (node_idx, i))
        .collect();

    // Same threshold and noise-tolerant span measurement as
    // find_structural_bubbles, so `BubbleSite.is_structural` always agrees
    // with the classification actually used for phasing -- see that
    // function's and materialize_arm_len_tolerant's doc comments for why.
    let threshold = ((n_reads as f64 * min_allele_freq).ceil() as i32).max(1);

    find_bubbles(nodes, topo, n_reads, min_allele_freq)
        .into_iter()
        .filter_map(|(entry_node, arm_starts)| {
            let consensus_pos = *path_pos.get(&entry_node)?;

            let arm_read_counts: Vec<u32> = arm_starts
                .iter()
                .map(|&a| {
                    edge_reads
                        .get(&(entry_node, a))
                        .map_or(0, |v| v.len() as u32)
                })
                .collect();

            // Mirror the guard in find_structural_bubbles: a 0-length arm
            // (arm_start == exit) is not structural regardless of what follows.
            let is_structural = arm_starts.iter().any(|&a| {
                let len = if real_in_edge_count(nodes, a, threshold) > 1 {
                    0
                } else {
                    materialize_arm_len_tolerant(nodes, a, phasing_bubble_min_span, threshold)
                };
                len >= phasing_bubble_min_span
            });

            let arm_sequences: Vec<Vec<u8>> =
                arm_starts.iter().map(|&a| arm_sequence(nodes, a)).collect();

            Some(BubbleSite {
                consensus_pos,
                arm_read_counts,
                arm_sequences,
                is_structural,
            })
        })
        .collect()
}

// ─── Structural-split validity check ──────────────────────────────────────────
//
// `find_structural_bubbles`'s own per-bubble vote-count threshold only checks
// whether a single bubble's arm split looks significant; it cannot tell a
// genuine second haplotype apart from a column that happens to split reads
// close to evenly by coincidence (a periodic-repeat phase-registration tie --
// Known Bugs #3/#4/#9's family). Confirmed on real GAA×30/GAA×100 data: once
// `find_structural_bubbles`'s noise-tolerant span measurement (above) started
// correctly finding real structural bubbles instead of none, `phasing_groups`
// began correctly refusing to bridge conflicting reads together (see its own
// doc comment) -- but that, in turn, revealed several small clusters driven by
// nothing more than a handful of reads agreeing on one coincidental tie
// column. A genuine second allele's read population should ALSO show a
// distinct read-length mode (the two alleles' own lengths differ) and hold up
// on its own vote weight; a coincidental tie column produces reads whose
// lengths are indistinguishable from the group they were spuriously split out
// of.

/// Minimum read-length gap (bp) between two candidate groups' medians before
/// a length-based split is trusted at all. Guards against a near-zero
/// absolute gap being amplified into an apparently "significant" separation
/// ratio purely because both groups happen to have very tight internal
/// spread (e.g. two 1-read groups always have zero spread).
const MIN_LENGTH_GAP_BP: f64 = 8.0;

/// A length gap must be at least this many multiples of the pooled spread
/// (median absolute deviation) to count as genuine bimodal separation rather
/// than ordinary scatter around one shared mode. 3 MADs is a conservative
/// bar (roughly analogous to ~2 standard deviations for a normal-ish
/// distribution): real independent haplotype length populations in the
/// validated CAG/GAA scenarios differ by tens to hundreds of bp against a
/// MAD on the order of single-digit bp (ordinary indel/substitution read
/// noise), while the confirmed spurious clusters' medians sit within a few
/// bp of the group they were actually noise from.
const LENGTH_SEPARATION_MADS: f64 = 3.0;

/// Floor on the per-group spread estimate, so two groups that both happen to
/// have very tight (near-identical) internal lengths -- including the
/// degenerate 1-read-group case, spread 0 -- don't produce an artificially
/// inflated separation ratio from a near-zero denominator.
const MIN_SPREAD_FLOOR_BP: f64 = 3.0;

/// Median of `vals`, treated as `f64`. Empty input returns 0.0 (callers only
/// ever pass non-empty read-length lists here).
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

/// Median absolute deviation of `vals` around `median`.
fn mad_usize(vals: &[usize], median: f64) -> f64 {
    let dev: Vec<usize> = vals
        .iter()
        .map(|&x| (x as f64 - median).abs().round() as usize)
        .collect();
    median_usize(&dev)
}

/// True when two groups' read-length distributions are genuinely bimodally
/// separated (real length-differing alleles) rather than ordinary scatter
/// around one shared mode (see the constants above for the exact bar).
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

/// Validates candidate haplotype groups (as produced by [`phasing_groups`],
/// before any of its own below-`min_reads` folding is undone) against a
/// minimum-depth vote-weight floor, plus -- when 2+ structural bubbles were
/// involved -- read-length bimodality as the primary discriminating signal.
/// Any candidate that fails merges into whichever already-accepted group its
/// own read lengths are closest to.
///
/// # Why the vote-weight bar is `min_reads`, not `min_allele_freq`
///
/// An earlier version used the same `min_allele_freq`-derived threshold
/// `find_bubbles`/`find_structural_bubbles` use to decide whether a single
/// bubble's own arm is significant. That is the wrong bar to apply to a
/// *cluster*'s final size here: `phasing_groups`'s whole point is to keep
/// only reads that agree *consistently across every bubble they have a
/// recorded choice at* (see its own doc comment on why plain union-find
/// bridging is unsound), which is a strictly harder bar to clear than any
/// single bubble's own arm weight -- a read that is noisy at just one of
/// several bubbles is correctly excluded from the clean cluster even though
/// it still supports the same haplotype overall. Reusing the un-attrited
/// `min_allele_freq * n_reads` threshold against this already-attrited count
/// double-penalises exactly the reads `phasing_groups` is supposed to filter
/// once, rejecting genuine minority alleles outright: confirmed on real data
/// (`multi_skew_cag20_40`, a true ~20%-frequency CAG×40 minority allele) where
/// the cleanly cross-bubble-consistent cluster (5 reads) cleared its own
/// bubbles' individual arm-weight thresholds fine but fell under the
/// pool-wide `min_allele_freq` bar purely from that attrition, even though it
/// was clearly read-length-separated from the majority allele. `min_reads`
/// (the same absolute depth floor `phasing_groups` already uses for ITS OWN
/// below-floor folding, and the same floor `PoaGraph` requires before
/// attempting a consensus at all) is the bar actually being asked here --
/// "is there enough depth to build a consensus from this cluster at all" --
/// so reusing it, rather than the separate significance-of-a-single-arm bar,
/// is the correct comparison.
///
/// # Why length separation is required only when 2+ bubbles are in play
///
/// A genuine haplotype split is frequently the SAME length on both arms by
/// construction -- not just same-length SNP-bubble-fallback haplotypes
/// (`structural_bubble_phasing_ignores_snp_bubbles`, handled by a completely
/// separate code path anyway), but same-length *structural* bubbles too: a
/// single long homopolymer run that differs by base rather than length
/// (confirmed by `locked_arm_deep_bubble_alleles_lost`, a G×30-vs-A×30 arm,
/// which regressed under an earlier version of this function that required
/// length separation unconditionally). Requiring length separation whenever
/// there is only one structural bubble in the whole graph would reject that
/// class of real call outright.
///
/// The confirmed failure mode this function targets -- a coincidental
/// phase-registration tie in a periodic repeat producing a column whose vote
/// split looks significant but carries no real length signal -- was only
/// ever observed (and is only plausible) when 2+ structural bubbles exist:
/// `find_structural_bubbles` evaluates every fork in the graph and keeps
/// only those whose (now noise-tolerant) span clears the threshold, so
/// multiple surviving bubbles means multiple independent candidates were
/// compared, and a coincidental near-even split becomes far more likely to
/// slip through *some* one of them than when there is only a single
/// candidate in the entire graph. Gating the length requirement on bubble
/// count is therefore not an arbitrary threshold but a direct response to
/// that multiple-candidates dynamic: one bubble carries no such selection
/// risk (its own span survived on genuine structural grounds, not on being
/// picked as the most-balanced-looking of many), so the depth floor alone is
/// trusted; two or more reintroduces exactly the risk length separation is
/// meant to catch.
///
/// # Why a clean distinguishing bubble is *also* required (Phase 3 of
/// `design/bypass_edge_delete_rework.md`)
///
/// Length separation alone became insufficient once Delete ops were
/// represented as bypass edges (pure bypass). A single allele's deletion-heavy
/// reads no longer register at the structural bubbles they skipped (their
/// bypass edges are invisible to `edge_reads`, which the per-read arm signature
/// is built from), and a deletion elsewhere frame-shifts a read onto the
/// *wrong* (shorter) arm at a downstream bubble. `phasing_groups` then splits
/// that one allele into sub-groups that happen to be narrowly length-separated
/// from each other -- confirmed on `multi_skew_cag20_40`, where the true CAG×40
/// allele over-split into ~39- and ~32-unit sub-groups whose 16bp median gap
/// just cleared the 12bp length-separation bar. Length separation cannot tell a
/// genuinely bimodal two-allele length distribution from a unimodal one that
/// `phasing_groups` bisected, because it only ever compares the two
/// already-tight halves.
///
/// The added requirement is orthogonal to length: a candidate is a genuinely
/// distinct allele only if it also has a **clean distinguishing bubble** versus
/// every already-confirmed group -- a structural bubble where a clear majority
/// (>= `CLEAR_MAJORITY`) of each group, over enough reads (>= `MIN_ARM_COV`),
/// takes a *different* arm. Measured (see the design doc's Phase 3 spike):
/// every genuine two-allele case -- `multi_gaa30_100`, `multi_cag15_25`,
/// `multi_cag20_50`, `structural_phasing_no_contamination_on_noisy_periodic_repeat`,
/// `structural_phasing_small_gap_bridge_candidate_stress` -- has such a bubble
/// at majority fractions 0.80-1.00, while `multi_skew_cag20_40`'s two CAG×40
/// sub-groups take the *same* majority arm at every bubble (no distinguishing
/// bubble at all), so they fail and merge back to one allele. A candidate that
/// clears length separation but fails this check is the same allele split by
/// deletion noise; it merges wholesale into the length-nearest confirmed group
/// it is structurally indistinguishable from. This check, like length
/// separation, is gated on 2+ structural bubbles (with 0-1 there is nothing to
/// distinguish on beyond length, and the SNP-bubble fallback is a separate
/// path). The `else` per-read length-routing branch (the `multi_gaa30_100`
/// contamination fix, for candidates that fail *length* separation) is
/// untouched.
fn validate_and_merge_groups(
    groups: Vec<Vec<usize>>,
    reads: &[Vec<u8>],
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

    // Clear-majority threshold and minimum arm coverage for the "clean
    // distinguishing bubble" check (validated in the Phase 3 spike; the
    // outcome is robust across 0.60-0.70, 0.60 is what was measured).
    const CLEAR_MAJORITY: f64 = 0.60;
    const MIN_ARM_COV: usize = 2;

    // Per-bubble `read -> arm index` map -- the same arm signal `phasing_groups`
    // clustered on, rebuilt here from `edge_reads` (cheap: O(reads x bubbles),
    // and only when 2+ bubbles gate the check below).
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

    // (majority arm, majority fraction, coverage) for a group at bubble `b`.
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

    // True iff some structural bubble has both groups' clear majorities on
    // *different* arms (with enough coverage on each side).
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

    // The largest group is always trusted as the baseline/reference -- same
    // philosophy as every other bubble-selection heuristic in this file
    // (e.g. partition_reads_by_bubble's own largest-group fallback target).
    let mut confirmed: Vec<Vec<usize>> = vec![groups[0].clone()];
    let mut confirmed_lens: Vec<Vec<usize>> =
        vec![groups[0].iter().map(|&r| reads[r].len()).collect()];

    for cand in groups.into_iter().skip(1) {
        let cand_lens: Vec<usize> = cand.iter().map(|&r| reads[r].len()).collect();

        let significant = cand.len() >= min_reads;
        let separated_from_all = !require_length_separation
            || confirmed_lens
                .iter()
                .all(|existing| length_separated(&cand_lens, existing));

        if significant && separated_from_all {
            // Length-separated -- but a genuinely distinct allele must ALSO
            // have a clean distinguishing bubble vs every confirmed group
            // (gated on 2+ bubbles, same as length separation). A candidate
            // that clears length but fails this is one allele's
            // deletion-noise-split sub-group; merge it wholesale into the
            // length-nearest confirmed group it can't be distinguished from.
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
            // Only one confirmed group exists so far -- nothing to route
            // individual members toward yet, merge wholesale (this is also
            // the common, harmless case: a small or non-separated candidate
            // whose members really do all belong together).
            confirmed_lens[0].extend(cand_lens);
            confirmed[0].extend(cand);
        } else {
            // Not trustworthy as its own allele. Route each *member*
            // individually to whichever confirmed group its own length is
            // closest to, rather than moving the whole candidate as one
            // block to wherever the candidate's *aggregate* median lands.
            //
            // Why per-read, not per-group: a rejected candidate can itself
            // be an internally mixed cluster -- confirmed on real
            // GAA×30/GAA×100 HiFi data (`multi_gaa30_100`): a 5-read
            // candidate contained 4 reads of one true allele and 1 read of
            // the other, all sharing the exact same cross-bubble arm-choice
            // pattern by coincidence (see `phasing_groups`'s doc comment on
            // this same underlying alignment-tie mechanism). The
            // candidate's own aggregate median length was pulled toward the
            // 4-read majority and the whole block -- including the 1
            // genuinely different read -- got merged into that majority's
            // confirmed group, misassigning the minority read. Each read's
            // own length is not similarly diluted, and the two confirmed
            // groups here are, by construction, already separated on
            // length (that is what qualified them as confirmed in the
            // first place), so per-read nearest-length routing recovers
            // the minority read correctly without affecting the majority.
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
                confirmed_lens[target].push(reads[r].len());
                confirmed[target].push(r);
            }
        }
    }

    confirmed
}

/// Groups reads by arm-choice compatibility across all structural bubbles.
///
/// For each structural bubble, edge_reads tells us which reads entered each arm.
/// Two reads are in the same haplotype group when they agree on every bubble
/// where both have a recorded arm choice. Reads with no arm assignment at any
/// structural bubble (pre-dating all bubbles, or not spanning them) are folded
/// into the largest group, as are reads whose partial information can't
/// distinguish between two or more existing groups (see below). Groups below
/// min_reads are also folded into the largest.
///
/// # Why this is not a plain pairwise union-find
///
/// An earlier version unioned any two reads whenever they never *conflicted*
/// (treating a missing arm choice at either read as an automatic pass) and
/// took the transitive closure of that relation. That is unsound whenever a
/// read has partial information (an arm choice at some bubbles but not
/// others): such a read can be pairwise-"compatible" with two *mutually
/// contradictory* fully-specified reads at once (e.g. read X entered arm 0 at
/// bubble A and has no recorded choice at bubble B; reads Y and Z both entered
/// arm 0 at bubble A but *disagree* at bubble B) -- plain union-find, being
/// purely pairwise, transitively merges Y and Z into one group through X even
/// though Y and Z themselves directly conflict and would never be unioned if
/// compared to each other. Confirmed on real data (a periodic GAA repeat
/// where `find_structural_bubbles`' noise-tolerant span measurement now
/// correctly detects 2-3 real structural bubbles instead of 0): reads with
/// only partial bubble coverage bridged two haplotypes with genuinely
/// conflicting votes into a single collapsed group, when they should have
/// stayed apart.
///
/// Fixed by building clusters incrementally, most-informative-read first,
/// and only ever merging a read into a cluster when it is compatible with
/// *exactly one* existing cluster -- never letting a read's own missing
/// information be used to bridge two clusters that would conflict with each
/// other. A read compatible with zero clusters starts a new one (this is
/// what correctly produces 3+ groups for nested bubbles, e.g. a short/
/// medium/long three-allele split: the short reads never reach the inner
/// bubble at all, so they are informationally *incompatible* with both the
/// medium and long clusters at the outer bubble -- a real conflict, not
/// ambiguity -- and correctly form their own third group). A read compatible
/// with two or more clusters carries no information that distinguishes them,
/// so -- exactly like a read with no recorded arm choice anywhere -- it is
/// folded into the largest final group rather than risking a wrong merge.
fn phasing_groups(
    edge_reads: &HashMap<(usize, usize), Vec<u32>>,
    bubbles: &[(usize, Vec<usize>)],
    n_reads: usize,
    min_reads: usize,
    reads: &[Vec<u8>],
) -> Vec<Vec<usize>> {
    let n_bubbles = bubbles.len();

    // sig[read][bubble] = Some(arm_idx) if that read entered that bubble arm, else None.
    let mut sig: Vec<Vec<Option<usize>>> = vec![vec![None; n_bubbles]; n_reads];
    for (b, (entry, arm_starts)) in bubbles.iter().enumerate() {
        for (arm_idx, &arm_start) in arm_starts.iter().enumerate() {
            if let Some(reads) = edge_reads.get(&(*entry, arm_start)) {
                for &r in reads {
                    let r = r as usize;
                    if r < n_reads {
                        sig[r][b] = Some(arm_idx);
                    }
                }
            }
        }
    }

    // Split reads into those with at least one recorded arm choice and those with none.
    let mut assigned: Vec<usize> = Vec::new();
    let mut unassigned: Vec<usize> = Vec::new();
    for (r, row) in sig.iter().enumerate() {
        if row.iter().any(|s| s.is_some()) {
            assigned.push(r);
        } else {
            unassigned.push(r);
        }
    }

    // Process the most-informative reads (fewest missing bubble choices)
    // first, ties broken by read index for determinism, so cluster
    // "templates" are established from the strongest evidence before any
    // partial-information read is considered.
    let mut order = assigned.clone();
    order.sort_by_key(|&r| {
        let known = sig[r].iter().filter(|s| s.is_some()).count();
        (std::cmp::Reverse(known), r)
    });

    // clusters[c][b] = this cluster's observed arm choice at bubble b (None
    // if no member has reached it yet). Grown incrementally: once a cluster
    // observes a value at a bubble, that becomes a hard requirement for
    // every subsequent read considered for that cluster.
    let mut clusters: Vec<Vec<Option<usize>>> = Vec::new();
    let mut cluster_members: Vec<Vec<usize>> = Vec::new();
    // Reads compatible with 2+ clusters at once: their partial information
    // can't tell those clusters apart on its own. Each carries the list of
    // clusters it's actually compatible with, so the resolution pass below
    // can pick among *those* rather than guessing globally.
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

    // Resolve bridge candidates and fully-unassigned reads against the
    // clusters formed above, using each read's own length as a
    // corroborating signal -- but only when that read's length is itself
    // plausible evidence, not an artifact of a partial (non-spanning) read.
    //
    // Why length at all: the previous design dumped both kinds of
    // leftover read unconditionally into whichever cluster was globally
    // largest, on the theory that "no unique bubble-arm signal" meant "no
    // signal at all," so the majority was the least-bad guess. That holds
    // for a read whose bubble-derived compatible set spans clusters that
    // all end up on the *same* true side. Confirmed wrong on real
    // GAA×30/GAA×100 HiFi data (`multi_gaa30_100`): a bridge candidate can
    // be compatible with clusters on *opposite* true sides purely because
    // one of its bubbles happens to overlap both (e.g. compatible with a
    // small cluster representing a genuine cross-bubble alignment tie --
    // see `validate_and_merge_groups`'s doc comment -- and also with the
    // far larger pure-majority cluster on the *other* allele). Dumping such
    // a read into the globally largest cluster is then actively wrong, not
    // a low-information default. Read length is a signal largely
    // orthogonal to bubble topology (a structural bubble is detected
    // *because* of a length difference), so it can break this kind of tie
    // even when the bubble evidence itself was misleading for this one
    // read.
    //
    // Why length is gated on plausibility, not applied unconditionally: a
    // first version of this fix compared every leftover read's length
    // against *all* clusters unconditionally and regressed immediately --
    // partial (non-spanning) reads that never reached any bubble have a
    // truncated length reflecting how much of the flank+repeat they
    // happened to cover, not their true allele's length, so several
    // partial reads of the *same* true allele clustered together purely
    // because they were all similarly short, fabricating a spurious extra
    // "allele" that `validate_and_merge_groups` then couldn't tell apart
    // from a real one (confirmed: `multi_gaa30_100` went from 2 output
    // alleles to a spurious 3rd, a 4-read all-partial, all-same-true-origin
    // group, once this fix's first version routed them by raw length
    // alone). A read's length is only trustworthy as a discriminating
    // signal when it is at least roughly as long as the *shortest*
    // well-supported (>= min_reads) cluster's own median -- i.e. plausible
    // as a genuine full instance of *some* known allele, not obviously
    // truncated relative to every one of them. `PLAUSIBLE_LEN_FRACTION`
    // (0.85) leaves room for a moderate deletion/error in an otherwise
    // full-length short-allele read without being so loose it accepts an
    // actually-partial read.
    //   - A bridge candidate that fails the plausibility check falls back
    //     to its own largest *bubble-compatible* cluster (still better than
    //     the global default: at least respects the partial bubble
    //     evidence it does have, which -- for a short/partial read -- is
    //     more trustworthy than its truncated length).
    //   - A fully unassigned read that fails the check has no compatible
    //     set to fall back to at all (zero recorded arm choices anywhere),
    //     so it defaults to the single largest cluster overall, exactly as
    //     before this fix.
    const PLAUSIBLE_LEN_FRACTION: f64 = 0.85;

    if !cluster_members.is_empty() {
        let snapshot_lens: Vec<f64> = cluster_members
            .iter()
            .map(|members| {
                median_usize(&members.iter().map(|&r| reads[r].len()).collect::<Vec<_>>())
            })
            .collect();
        let snapshot_sizes: Vec<usize> = cluster_members.iter().map(|m| m.len()).collect();
        let min_full_len = (0..cluster_members.len())
            .filter(|&c| snapshot_sizes[c] >= min_reads)
            .map(|c| snapshot_lens[c])
            .fold(f64::INFINITY, f64::min);

        let nearest_cluster_among = |read_idx: usize, candidates: &[usize]| -> usize {
            let len = reads[read_idx].len() as f64;
            candidates
                .iter()
                .copied()
                .min_by(|&a, &b| {
                    let da = (snapshot_lens[a] - len).abs();
                    let db = (snapshot_lens[b] - len).abs();
                    da.partial_cmp(&db)
                        .unwrap()
                        // Tie-break toward the larger (more-evidenced) cluster, then
                        // lowest index, for determinism.
                        .then_with(|| snapshot_sizes[b].cmp(&snapshot_sizes[a]))
                        .then_with(|| a.cmp(&b))
                })
                .unwrap()
        };
        let largest_overall = (0..cluster_members.len())
            .max_by_key(|&c| snapshot_sizes[c])
            .unwrap();
        let all_clusters: Vec<usize> = (0..cluster_members.len()).collect();
        // Fallback for a bridge candidate whose own length isn't trustworthy:
        // pick the largest cluster among its bubble-compatible set *by
        // member count*, not by length -- a read whose length we've just
        // decided not to trust must not turn around and use that same
        // length to break the tie among its fallback candidates.
        let largest_compatible = |candidates: &[usize]| -> usize {
            *candidates
                .iter()
                .max_by_key(|&&c| snapshot_sizes[c])
                .unwrap()
        };

        let mut resolved: Vec<(usize, usize)> = Vec::new();
        for (r, compat) in &bridge_candidates {
            let len = reads[*r].len() as f64;
            let target = if min_full_len.is_finite() && len >= PLAUSIBLE_LEN_FRACTION * min_full_len
            {
                nearest_cluster_among(*r, &all_clusters)
            } else {
                largest_compatible(compat)
            };
            resolved.push((*r, target));
        }
        for &r in &unassigned {
            let len = reads[r].len() as f64;
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
        // No cluster exists at all (no read anywhere had a recorded arm
        // choice) -- nothing to compare lengths against, so every read is
        // definitionally in one, undifferentiated group.
        cluster_members.push(Vec::new());
        cluster_members[0].extend(bridge_candidates.into_iter().map(|(r, _)| r));
        cluster_members[0].extend(unassigned);
    }

    let mut groups: Vec<Vec<usize>> = cluster_members;
    groups.sort_unstable_by_key(|g| std::cmp::Reverse(g.len()));

    if groups.is_empty() {
        groups.push(Vec::new());
    }

    // Fold groups below min_reads into the largest group.
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
            // Seed's own bases are genuine confirmation, not a skip-through.
            edge_reads.insert((i, i + 1), vec![0u32]);
        }

        Ok(PoaGraph {
            nodes,
            config,
            n_reads: 1,
            reads: vec![seed.to_vec()],
            edge_reads,
            edge_delete_reads: HashMap::new(),
            bypass_edges: HashMap::new(),
            warnings: 0,
            cached_spine: Vec::new(),
            spine_updated_at: 0,
            spine_interval: 1,
            align_scratch: AlignScratch::new(),
            spine_mers: HashMap::new(),
            fork_arm_index: HashMap::new(),
            used_band_retry: false,
        })
    }

    pub fn add_read(&mut self, read: &[u8]) -> Result<(), PoaError> {
        if read.is_empty() {
            return Err(PoaError::EmptyInput);
        }

        let (topo, rank_of) = topological_order(&self.nodes);
        // Permanent internal invariant, zero-cost in release builds: the
        // graph must always be a DAG entering `add_read`. If it isn't,
        // `topological_order`'s Kahn's-algorithm queue silently excludes the
        // cyclic nodes from `topo` rather than erroring, which would hang
        // `heaviest_path`'s predecessor-pointer traceback in an unbounded
        // loop the next time the spine is refreshed (see
        // `reuse_would_create_back_edge`'s doc comment for how this was
        // introduced and fixed). Catching the violation here, right at the
        // point a caller could next observe it, is far cheaper to diagnose
        // than the OOM it would otherwise cause.
        debug_assert_eq!(
            topo.len(),
            self.nodes.len(),
            "graph has a cycle before adding read {}: topo includes {} of {} nodes",
            self.n_reads,
            topo.len(),
            self.nodes.len()
        );

        // Refresh the spine when the cache is empty or the interval has elapsed.
        // The interval doubles each time the spine is stable (≤ SPINE_STABLE_THRESHOLD
        // base changes), and resets to 1 when it changes significantly.  The final
        // consensus always recomputes the heaviest path from scratch.  A stale
        // spine affects both alignment speed and anchor quality: shifted k-mers
        // built from an outdated spine can produce wrong anchor chains.
        let reads_since_update = self.n_reads.saturating_sub(self.spine_updated_at);
        if self.cached_spine.is_empty() || reads_since_update >= self.spine_interval {
            let new_spine = heaviest_path(&self.nodes, &topo, &rank_of, &self.bypass_edges);
            let diff = spine_diff(&self.cached_spine, &new_spine);
            self.cached_spine = new_spine;
            self.spine_updated_at = self.n_reads;
            // Rebuild spine minimizer index whenever the spine is refreshed.
            self.spine_mers = build_spine_mers(&self.cached_spine, MINI_K, MINI_W);
            if diff <= SPINE_STABLE_THRESHOLD {
                self.spine_interval = (self.spine_interval * 2).min(SPINE_MAX_INTERVAL);
            } else {
                self.spine_interval = 1;
            }
        }

        // Build minimizer anchor chain: match read k-mers against unique spine
        // k-mers, then retain the longest colinear hit chain via LIS.  The chain
        // constrains the per-row j-window inside align(), narrowing the DP band
        // on spine nodes at exact anchor positions.
        //
        // Density gate: in repetitive sequence (CAG, GAA, AAGGG repeats), error-
        // induced unique k-mers in the spine can match error k-mers in a read at
        // a different repeat unit but the same cycle phase, creating wrong anchors
        // with plausible (read_pos, topo_rank) pairs that survive LIS colinearity
        // filtering.  These wrong anchors narrow the j-window and may exclude the
        // correct j position, causing repeat unit loss.
        //
        // A fixed minimum threshold (MINI_MIN_CHAIN) separates the two regimes:
        // in repetitive sequence, expected wrong anchor count ≈ n_err² / (n_units
        // × n_sub_types) ≈ 1-3 for typical ONT depths — well below the threshold.
        // In non-repetitive sequence, ~46% of k-mers survive unmodified at 5%
        // error rate, giving ≈read_len × 0.46 / MINI_W correct anchors — well
        // above the threshold for reads ≥ 200 bp.
        let read_mers = compute_minimizers(read, MINI_K, MINI_W);
        let raw_anchors =
            build_anchor_chain(&read_mers, &self.spine_mers, &self.cached_spine, &rank_of);
        let anchors: Vec<(usize, usize)> = if raw_anchors.len() >= MINI_MIN_CHAIN {
            raw_anchors
        } else {
            vec![]
        };

        let (ops, retried) = align_with_retry(
            &self.nodes,
            &topo,
            &rank_of,
            &self.cached_spine,
            read,
            &self.config,
            &mut self.align_scratch,
            &anchors,
        )?;
        if retried {
            self.used_band_retry = true;
        }
        let read_idx = self.n_reads as u32;
        add_to_graph(
            &mut self.nodes,
            &mut self.edge_reads,
            &mut self.bypass_edges,
            &mut self.fork_arm_index,
            &rank_of,
            read,
            &ops,
            read_idx,
        );
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
            let mut graph_stats =
                compute_stats(&self.nodes, self.config.min_allele_freq, self.n_reads);
            graph_stats.median_input_read_len = median_read_len(&self.reads);
            return Ok(Consensus {
                sequence,
                coverage,
                path_weights,
                n_reads: 1,
                graph_stats,
                gaps: vec![],
                bubble_sites: vec![],
                read_indices: vec![],
            });
        }

        let (topo, rank_of) = topological_order(&self.nodes);

        let filtered: Vec<(usize, u8, i32)> = match self.config.consensus_mode {
            ConsensusMode::HeaviestPath => {
                let path = heaviest_path(&self.nodes, &topo, &rank_of, &self.bypass_edges);

                // Boundary trim: find the first/last node whose Match
                // coverage clears an absolute population floor. This has
                // to stay absolute, not evidence-relative -- a trailing (or
                // leading) minority extension has delete_count == 0 at its
                // first node (nobody deletes through it; the rest of the
                // population simply never reaches that far under
                // semi-global alignment), so a Match-vs-Delete comparison
                // alone can't distinguish "real interior disagreement"
                // from "the majority already ended here." Only an absolute
                // floor can.
                //
                // The floor's *population* is a per-position local estimate
                // (see `local_population_profile`), not the flat
                // `self.n_reads` used before: for a genuinely partial-read
                // population (no single read spans the whole target
                // region), `self.n_reads` counts reads that can never
                // reach a given position at all, making
                // `(n_reads/2+1).max(2)` unreachable almost everywhere and
                // silently collapsing the consensus toward half its true
                // length (confirmed on `gen_short_reads_long_region_ont_r10`
                // in `bench/compare_callers.py --general`).
                //
                // Radius choice: see `LOCAL_POP_RADIUS`'s own doc comment
                // for the full empirical trade-off (a read-length-scaled
                // radius under-fixes this bug; a small, jitter-scale radius
                // over-fixes it by resurrecting three majority-vs-minority
                // regressions on real repeat data). For a normal
                // fully-spanning population this is a no-op (nearby columns
                // already share the same near-total coverage, same as
                // before); only a structurally partial population sees the
                // floor actually shrink.
                let radius = LOCAL_POP_RADIUS;
                let path_coverage: Vec<u32> = path
                    .iter()
                    .map(|&(node_idx, _, _)| self.nodes[node_idx].coverage)
                    .collect();
                let local_pop = local_population_profile(&path_coverage, radius);
                let local_min_cov_by_pos: Vec<u32> = local_pop
                    .iter()
                    .map(|&pop| coverage_threshold(pop as usize, self.config.min_coverage_fraction))
                    .collect();

                let meets_floor: Vec<bool> = path
                    .iter()
                    .enumerate()
                    .map(|(i, &(node_idx, _, _))| {
                        self.nodes[node_idx].coverage >= local_min_cov_by_pos[i]
                    })
                    .collect();

                let start = meets_floor.iter().position(|&s| s).unwrap_or(0);
                let end = meets_floor
                    .iter()
                    .rposition(|&s| s)
                    .map(|i| i + 1)
                    .unwrap_or(path.len());

                let range = if start < end {
                    start..end
                } else {
                    0..path.len()
                };

                // Interior inclusion: within the trimmed boundary, two
                // distinct evidence axes decide whether a node survives,
                // and they need different tests because they answer
                // different questions:
                //
                // 1. Match vs Delete -- "is this base here at all, or did
                //    the reads that reached this exact point skip it?"
                //    Only meaningful when NO fork (2+ out-edges) exists
                //    anywhere back along this node's unbranched run, all
                //    the way to the trimmed boundary: with no competing
                //    arm anywhere in reach, a node's own coverage and
                //    delete_count are the *exact*, complete count of reads
                //    that reached this precise point, so a plain
                //    `coverage > delete_count` is sufficient and needs no
                //    population floor -- majority-delete fabrication and a
                //    false rescue of it (Known Bugs #6, #9; DAB1 SCA37
                //    cov=3 vs del=23) both come down to exactly this axis.
                // 2. This base vs a *different* base -- a real fork exists
                //    somewhere back along this node's unbranched run. Here
                //    "coverage beats delete_count" is the wrong test: a
                //    base can have zero deletes of its own and still be
                //    nothing more than one arm of a genuine near-tie with
                //    a sibling arm (confirmed regressions: RFC1 AAAAG
                //    5-vs-5, SV tandem-duplication 10-vs-4 -- the latter
                //    only visible several nodes downstream of the fork
                //    itself, which is why the search below walks the whole
                //    unbranched run rather than checking one hop back).
                //    This axis needs the fork's own coverage to clear a
                //    majority of the fork's local total (not merely beat
                //    zero), gated on that local total itself clearing the
                //    position-local absolute floor (see
                //    `local_population_profile`) so a heavily fragmented,
                //    near-empty fork can't manufacture a "majority" out of
                //    noise (Known Bug #7). This gate used a single global
                //    `min_cov` before the partial-read-population fix
                //    below; using the same per-position local floor here
                //    keeps both checks on one consistent notion of
                //    "population," and is a no-op for the fully-spanning
                //    populations Bug #7 was originally about (the nearby
                //    window still sees the same near-total coverage).
                //
                // A node several hops downstream of a fork, on an
                // otherwise-unbranched continuation of one of its arms,
                // still belongs to axis 2 -- the fork's accept/reject
                // verdict has to hold for its entire arm, not just the
                // arm's first node. Bug #8 originally found this via a
                // bounded (64-hop) backward walk over `path[]` re-derived on
                // every call; Phase 4 (design/graph_data_model_rework.md)
                // replaces that re-derivation with `Node.nearest_fork`, a
                // cache populated incrementally as the graph is built (see
                // `propagate_fork_if_new`) instead of walked fresh here.
                // Judgment is unchanged -- same two axes, same tests below --
                // only how the fork-context lookup happens.
                //
                // One boundary-compatibility check is still needed: the old
                // walk refused to look further back than `range_start` (the
                // start of *this call's* boundary-trimmed region), so a
                // fork sitting entirely within the trimmed-away leading
                // section was treated as "none found," falling back to
                // axis 1. `nearest_fork` is a graph-level cache with no
                // notion of any particular call's boundary trim, so that
                // check is reproduced explicitly via `rank_of` to keep the
                // exact same behaviour at that edge.
                let range_start = range.start;
                let range_start_rank = rank_of[path[range_start].0];
                range
                    .filter(|&i| {
                        let (node_idx, _, _) = path[i];
                        if meets_floor[i] {
                            return true;
                        }
                        let fork_info = self.nodes[node_idx]
                            .nearest_fork
                            .filter(|&(pred_idx, _)| rank_of[pred_idx] >= range_start_rank);
                        let Some((pred_idx, arm_idx)) = fork_info else {
                            return self.nodes[node_idx].coverage
                                > self.nodes[node_idx].delete_count;
                        };
                        // Phase 4 territory (design/graph_data_model_rework.md):
                        // this is the interior filter's fork-population/
                        // plurality logic (Known Bugs #6-#10). Kept on total
                        // (matched + deleted) weight, byte-identical to
                        // pre-split behaviour -- not part of this pass.
                        let local_total: i32 = self.nodes[pred_idx]
                            .out_edges
                            .iter()
                            .map(|&(_, ew)| ew.total())
                            .sum();
                        if (local_total.max(0) as u32) < local_min_cov_by_pos[i] {
                            return false;
                        }
                        let local_min_cov = coverage_threshold(
                            local_total.max(0) as usize,
                            self.config.min_coverage_fraction,
                        );
                        if self.nodes[node_idx].coverage as i32 >= local_min_cov as i32 {
                            return true;
                        }
                        // Plurality relaxation: trust the arm's own entry
                        // weight against its siblings' -- but *only* when
                        // both the fork itself and this candidate are
                        // "clean" (zero delete_count), i.e. genuinely
                        // undisputed populations on both ends, not entangled
                        // with an unresolved Match-vs-Delete decision. A
                        // fork with its own delete_count > 0 means the
                        // arrival at the fork is itself still contested
                        // (confirmed case: RFC1 AAAAG's fork had cov=4,
                        // del=6 of its own -- a majority-delete decision one
                        // level up from the 5-vs-5 split below it -- so
                        // trusting that split's plurality independently of
                        // the unresolved arrival fabricated a base no read
                        // actually has). Confirmed target case: a genuine
                        // 3-way split of the *full* population at a fork
                        // with cov=7, del=0 of its own (cag20_d05_r10) --
                        // no attrition, no unresolved arrival, just three
                        // different bases competing at the same position.
                        if self.nodes[pred_idx].delete_count != 0
                            || self.nodes[node_idx].delete_count != 0
                        {
                            return false;
                        }
                        let arm_weight = self.nodes[pred_idx]
                            .out_edges
                            .iter()
                            .find(|&&(to, _)| to == arm_idx)
                            .map(|&(_, ew)| ew.total())
                            .unwrap_or(0);
                        self.nodes[pred_idx]
                            .out_edges
                            .iter()
                            .all(|&(_, ew)| ew.total() <= arm_weight)
                    })
                    .map(|i| path[i])
                    .collect()
            }
            // `MajorityFrequency` keeps the flat, whole-population floor:
            // its documented use case is HiFi data with near-identical read
            // lengths within an allele group (CLAUDE.md's "Consensus
            // Modes"), not the genuinely-partial-read populations the local
            // profile above targets, and `topo` order (unlike `path`) does
            // not correspond to sequence position closely enough for a
            // position-windowed estimate to be meaningful.
            ConsensusMode::MajorityFrequency => {
                majority_frequency(&self.nodes, &topo, self.min_coverage())
            }
        };

        let sequence: Vec<u8> = filtered.iter().map(|&(_, base, _)| base).collect();
        let coverage: Vec<u32> = filtered
            .iter()
            .map(|&(node_idx, _, _)| self.nodes[node_idx].coverage)
            .collect();
        let path_weights: Vec<i32> = filtered.iter().map(|&(_, _, w)| w).collect();

        let mut graph_stats = compute_stats(&self.nodes, self.config.min_allele_freq, self.n_reads);
        graph_stats.median_input_read_len = median_read_len(&self.reads);
        let gaps = detect_coverage_gaps(&coverage);
        let bubble_sites = collect_bubble_sites(
            &self.nodes,
            &topo,
            &filtered,
            &self.edge_reads,
            self.config.min_allele_freq,
            self.n_reads,
            self.config.phasing_bubble_min_span,
        );

        Ok(Consensus {
            sequence,
            coverage,
            path_weights,
            n_reads: self.n_reads,
            graph_stats,
            gaps,
            bubble_sites,
            read_indices: vec![],
        })
    }

    pub fn stats(&self) -> GraphStats {
        let mut s = compute_stats(&self.nodes, self.config.min_allele_freq, self.n_reads);
        s.median_input_read_len = median_read_len(&self.reads);
        s
    }

    /// Number of long-unbanded warnings emitted during `add_read` calls.
    /// With the new bubble-aware aligner this always returns 0 (warnings removed).
    pub fn warnings_emitted(&self) -> usize {
        self.warnings
    }

    /// Return all edge weights in the graph (one entry per directed edge).
    ///
    /// Total (matched + deleted) traversal count, unchanged in meaning from
    /// before the Match/Delete edge-weight split -- this is a general
    /// diagnostic/visualization API (see `src/plot.rs`), not a consensus
    /// decision point, so its contract is preserved as-is.
    pub fn edge_weights(&self) -> Vec<i32> {
        self.nodes
            .iter()
            .flat_map(|n| n.out_edges.iter().map(|&(_, ew)| ew.total()))
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
        let spine = heaviest_path(&self.nodes, &topo, &rank_of, &self.bypass_edges);
        let (ops, _retried) = align_with_retry(
            &self.nodes,
            &topo,
            &rank_of,
            &spine,
            read,
            &self.config,
            &mut AlignScratch::new(),
            &[],
        )?;
        Ok((ops, 0, rank_of))
    }

    /// Like [`align_read_ops`] but kept for compatibility.
    pub fn align_read_ops_unbanded(
        &self,
        read: &[u8],
    ) -> Result<(Vec<AlignOp>, Vec<usize>), PoaError> {
        let (topo, rank_of) = topological_order(&self.nodes);
        let spine = heaviest_path(&self.nodes, &topo, &rank_of, &self.bypass_edges);
        let (ops, _retried) = align_with_retry(
            &self.nodes,
            &topo,
            &rank_of,
            &spine,
            read,
            &self.config,
            &mut AlignScratch::new(),
            &[],
        )?;
        Ok((ops, rank_of))
    }

    /// Number of nodes currently in the graph.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// `true` if any `add_read` call so far needed `align_with_retry`'s
    /// pass 2 or 3 (the configured band was too narrow for at least one
    /// read). Crate-internal: consumed by `build_graph` (`src/lib.rs`) to
    /// decide whether the whole graph needs a consistent unbanded rebuild
    /// -- see `used_band_retry`'s own field doc comment for why a graph
    /// built from a mix of band widths across different reads is not
    /// trustworthy as-is in a periodic/repetitive locus.
    pub(crate) fn used_band_retry(&self) -> bool {
        self.used_band_retry
    }

    /// Return a snapshot of the graph topology for visualization and inspection.
    ///
    /// Nodes are in topological order; `edges` use topological ranks as source
    /// and target identifiers.  `spine_ranks` lists the ranks of nodes on the
    /// heaviest-path (consensus) spine.
    pub fn graph_topology(&self) -> crate::types::GraphTopology {
        use crate::types::{GraphEdgeInfo, GraphNodeInfo, GraphTopology};

        let (topo, rank_of) = topological_order(&self.nodes);
        let spine = heaviest_path(&self.nodes, &topo, &rank_of, &self.bypass_edges);

        let nodes: Vec<GraphNodeInfo> = topo
            .iter()
            .enumerate()
            .map(|(rank, &node_idx)| {
                let n = &self.nodes[node_idx];
                GraphNodeInfo {
                    node_idx,
                    base: n.base,
                    coverage: n.coverage,
                    delete_count: n.delete_count,
                    topo_rank: rank,
                }
            })
            .collect();

        // GraphEdgeInfo.weight is documented as "number of reads that
        // traversed this edge" -- total (matched + deleted), unchanged from
        // before the split.
        let edges: Vec<GraphEdgeInfo> = topo
            .iter()
            .enumerate()
            .flat_map(|(from_rank, &node_idx)| {
                self.nodes[node_idx]
                    .out_edges
                    .iter()
                    .map(move |&(succ_idx, ew)| (from_rank, succ_idx, ew.total()))
                    .collect::<Vec<_>>()
            })
            .map(|(from_rank, succ_idx, weight)| GraphEdgeInfo {
                from_rank,
                to_rank: rank_of[succ_idx],
                weight,
            })
            .collect();

        let spine_ranks: Vec<usize> = spine
            .iter()
            .map(|(node_idx, _, _)| rank_of[*node_idx])
            .collect();

        GraphTopology {
            nodes,
            edges,
            spine_ranks,
        }
    }

    /// For each node, return `(out_edges_count, max_out_edge_weight, min_out_edge_weight)`.
    ///
    /// Total (matched + deleted) weight, same reasoning as `edge_weights()`.
    pub fn node_out_edge_info(&self) -> Vec<(usize, i32, i32)> {
        self.nodes
            .iter()
            .map(|n| {
                let count = n.out_edges.len();
                let max_w = n
                    .out_edges
                    .iter()
                    .map(|&(_, ew)| ew.total())
                    .max()
                    .unwrap_or(0);
                let min_w = n
                    .out_edges
                    .iter()
                    .map(|&(_, ew)| ew.total())
                    .min()
                    .unwrap_or(0);
                (count, max_w, min_w)
            })
            .collect()
    }

    /// Return arm lengths for each bubble node with 2+ qualifying arms (weight >= threshold).
    ///
    /// Arm length is measured by walking the single-successor chain from each arm start,
    /// stopping at reconvergence points (nodes with multiple in-edges after the first step).
    ///
    /// `weight_threshold` is compared against total (matched + deleted) weight,
    /// unchanged from before the Match/Delete edge-weight split -- a general
    /// diagnostic API (see `tests/rfc1_real_data.rs`), not a consensus decision point.
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
                .filter(|&&(_, ew)| ew.total() >= weight_threshold)
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
    ///
    /// Structural bubbles (arm span ≥ phasing_bubble_min_span) are phased first using
    /// cross-bubble compatibility grouping. If no structural bubbles are found and the
    /// alignment was banded, the graph is rebuilt with unbanded alignment so that large
    /// allele-length differences create a visible Insert-arm bubble; structural phasing
    /// then retries on the unbanded graph. If still no structural bubble, falls back to
    /// single-best SNP bubble partitioning for substitution haplotypes.
    ///
    /// Each returned [`Consensus`]'s `read_indices` are in this graph's own
    /// read ordering: the [`PoaGraph::new`] seed is index 0 and each
    /// [`PoaGraph::add_read`] call follows in order. The free-function
    /// [`crate::consensus_multi`] wrapper translates these back to the input
    /// slice's indexing; when calling this method directly the add-order
    /// indexing is left as-is.
    pub fn consensus_multi(&self) -> Result<Vec<Consensus>, PoaError> {
        if self.n_reads < self.config.min_reads {
            return Err(PoaError::InsufficientDepth {
                got: self.n_reads,
                min: self.config.min_reads,
            });
        }

        let (topo, _) = topological_order(&self.nodes);

        // Try structural bubble phasing first (length variants, SVs, large indels).
        let structural = find_structural_bubbles(&self.nodes, &topo, self.n_reads, &self.config);

        // If no structural bubble was found and alignment was banded, the band may be too
        // narrow to create a clean Insert-arm bubble for large allele-length differences.
        // Rebuild the graph with unbanded alignment (reads include flanking → rotation is
        // anchored) and recurse once. The config's band settings are zeroed so the recursive
        // call does not loop.
        if structural.is_empty() && (self.config.band_width > 0 || self.config.adaptive_band) {
            let mut cfg2 = self.config.clone();
            cfg2.band_width = 0;
            cfg2.adaptive_band = false;
            let seed = &self.reads[0];
            let mut g2 = PoaGraph::new(seed, cfg2)?;
            for read in self.reads.iter().skip(1) {
                g2.add_read(read)?;
            }
            return g2.consensus_multi();
        }

        let groups = if !structural.is_empty() {
            let g = phasing_groups(
                &self.edge_reads,
                &structural,
                self.n_reads,
                self.config.min_reads,
                &self.reads,
            );
            // Validate against read-length bimodality + vote weight + a clean
            // distinguishing bubble before trusting phasing_groups' split --
            // see validate_and_merge_groups' doc comment for why this is scoped
            // to the structural-bubble path specifically, not the same-length
            // SNP-bubble fallback below.
            validate_and_merge_groups(
                g,
                &self.reads,
                self.config.min_reads,
                &structural,
                &self.edge_reads,
            )
        } else {
            // Unbanded alignment and still no structural bubble. Try SNP bubble
            // partitioning for substitution haplotypes (same-length alleles).
            let snp_bubbles = find_bubbles(
                &self.nodes,
                &topo,
                self.n_reads,
                self.config.min_allele_freq,
            );
            if !snp_bubbles.is_empty() {
                let (entry, arm_starts) = snp_bubbles
                    .iter()
                    .max_by_key(|(entry, arms)| {
                        arms.iter()
                            .filter_map(|&arm| self.edge_reads.get(&(*entry, arm)))
                            .map(|v| v.len())
                            .min()
                            .unwrap_or(0)
                    })
                    .unwrap();
                partition_reads_by_bubble(&self.edge_reads, *entry, arm_starts, self.n_reads)
            } else {
                return Ok(vec![self.consensus()?]);
            }
        };

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
            let mut c = sub.consensus()?;
            c.read_indices = group.to_vec();
            results.push(c);
        }

        Ok(results)
    }

    fn min_coverage(&self) -> u32 {
        coverage_threshold(self.n_reads, self.config.min_coverage_fraction)
    }
}

/// Minimum coverage required to keep a node, given a population size (either
/// the full read count, or the local vote total at one bubble) and the
/// configured `min_coverage_fraction`.
fn coverage_threshold(population: usize, min_coverage_fraction: f64) -> u32 {
    if min_coverage_fraction > 0.0 {
        ((population as f64 * min_coverage_fraction).ceil() as u32).max(1)
    } else if population <= 1 {
        1
    } else {
        ((population / 2 + 1).max(2)) as u32
    }
}

/// For each position `i` along `coverages`, the maximum value observed within
/// `radius` positions on either side (inclusive) -- an O(n) sliding-window
/// maximum (monotonic deque), safe for path lengths in the thousands.
///
/// This is the "local population" proxy consumed by `coverage_threshold`'s
/// absolute-floor fallback for `ConsensusMode::HeaviestPath` (boundary trim
/// and the interior filter's fork-population gate). It replaces `self.n_reads`
/// -- the total read count ever added to the graph -- with an estimate of
/// how many reads could plausibly reach *near* this exact position.
///
/// The distinction matters for a genuinely partial-read population (no
/// single read spans the whole target region, by construction -- see
/// `Scenario.partial` / `gen_short_reads_long_region_ont_r10` in
/// `bench/compare_callers.py`): `self.n_reads` counts reads that structurally
/// can never reach a given interior position at all, so the absolute floor
/// `(n_reads/2+1).max(2)` is unreachable almost everywhere and the consensus
/// silently collapses toward half its true length. A *local* window instead
/// asks "how much agreement has this graph achieved *anywhere nearby*" --
/// for a fully-spanning population that is still essentially `n_reads`
/// (nearby columns share the same near-total coverage), so the fallback is
/// unchanged in the common case; for a partial population it correctly
/// shrinks to the achievable local depth.
///
/// A *global* (whole-path) maximum was tried first and rejected: a single
/// anomalously dense stretch elsewhere in the graph (e.g. where many
/// independently-placed partial reads' random start/end positions happen to
/// overlap) would inflate the floor everywhere else in the same graph,
/// which is exactly the miscalibration being fixed, just moved to a
/// different spot. Bounding the window to `radius` (see `LOCAL_POP_RADIUS`'s
/// own doc comment for how that constant's value was chosen -- it is an
/// empirically-tuned trade-off, not derived from read length or any other
/// single principle) keeps the estimate local to what nearby evidence can
/// actually support.
///
/// Deliberately does *not* use a forward-recovery scan (does coverage climb
/// back up further along?) -- Known Bug #9's writeup found that approach
/// regressed genuine trailing-minority-extension tests
/// (`long_repeat_length_majority_wins`, `edge_extreme_length_variation_majority_wins`)
/// by letting a stray noisy read overlapping the tail look like "recovery."
/// A symmetric, bounded window has no such directional bias: at a true
/// boundary, the high-coverage core sits on one side only, so the window's
/// max there still reflects the real population, correctly rejecting a
/// trailing extension a few positions further out.
fn local_population_profile(coverages: &[u32], radius: usize) -> Vec<u32> {
    let n = coverages.len();
    let mut out = vec![0u32; n];
    if n == 0 {
        return out;
    }
    let mut deque: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    let mut right = 0usize;
    for (i, out_i) in out.iter_mut().enumerate() {
        let hi = (i + radius).min(n - 1);
        while right <= hi {
            while let Some(&back) = deque.back() {
                if coverages[back] <= coverages[right] {
                    deque.pop_back();
                } else {
                    break;
                }
            }
            deque.push_back(right);
            right += 1;
        }
        let lo = i.saturating_sub(radius);
        while let Some(&front) = deque.front() {
            if front < lo {
                deque.pop_front();
            } else {
                break;
            }
        }
        *out_i = deque.front().map(|&idx| coverages[idx]).unwrap_or(0);
    }
    out
}

// ─── White-box tests: nearest_fork cache (Phase 4) ──────────────────────────
//
// `Node.nearest_fork` is a private field, so exercising it directly (not just
// observing its downstream effect on `consensus()`'s output) requires a
// same-module test, unlike the rest of this crate's tests (in `src/tests.rs`,
// a sibling module that only sees public API). This module is specifically
// for verifying the cache's own internal correctness -- the staleness
// question the design doc calls out as unproven.
#[cfg(test)]
mod fork_cache_tests {
    use super::*;

    fn b(s: &str) -> Vec<u8> {
        s.as_bytes().to_vec()
    }

    /// Directly exercises the staleness scenario: several reads establish a
    /// clean run of single-out-edge nodes, *then* a later read introduces a
    /// divergence at a node that was already a stable single-successor
    /// predecessor for existing descendants, turning it into a fork after
    /// the fact. Confirms the *pre-existing* descendants' cached
    /// `nearest_fork` gets updated by forward propagation, by comparing
    /// against an independent, exhaustive (unbounded) backward walk computed
    /// fresh from the raw graph -- not just checking the cache agrees with
    /// itself.
    #[test]
    fn nearest_fork_updates_for_preexisting_descendants_after_late_fork() {
        // Seed: 4 A's (unique-enough prefix) + a 12-base non-repetitive tail.
        // Node at position 4 (the tail's first base) is the one we'll fork
        // later; positions 5..15 are its pre-existing descendants.
        let seed = b("AAAACGTACGTACGTA");
        let cfg = PoaConfig {
            min_reads: 3,
            band_width: 0,
            adaptive_band: false,
            warn_on_long_unbanded: false,
            ..Default::default()
        };
        let mut g = PoaGraph::new(&seed, cfg).unwrap();

        // Several reads that match the seed exactly: establishes positions
        // 5..15 as a clean, unforked single-successor chain (their
        // nearest_fork should all be None at this point -- no fork anywhere
        // in this graph yet).
        for _ in 0..4 {
            g.add_read(&seed).unwrap();
        }
        for idx in 5..seed.len() {
            assert_eq!(
                g.nodes[idx].nearest_fork, None,
                "node {idx} should have no fork ancestor yet (pre-divergence)"
            );
            assert_eq!(
                g.nodes[idx].in_edges.len(),
                1,
                "node {idx} should be a plain single-predecessor chain node"
            );
        }
        let fork_node = 4; // 'C', the first tail base
        assert_eq!(g.nodes[fork_node].out_edges.len(), 1, "not yet a fork");
        let orig_child = g.nodes[fork_node].out_edges[0].0;
        assert_eq!(orig_child, 5);

        // Now a later read diverges immediately after position 4 (matches
        // "AAAAC" then substitutes the next base), giving node 4 a second
        // out-edge -- turning it into a fork *after* positions 5..15 already
        // existed as node 4's only descendants. Deliberately ends right at
        // the mismatch (no shared suffix appended): the read's remaining
        // bases would otherwise happen to equal seed[6..] exactly (both are
        // the same repeat unit), so the aligner would legitimately re-Match
        // them onto nodes 6..15 and reconverge there for real -- a correct
        // but different scenario (fork-then-genuine-reconvergence) that this
        // test isn't targeting. Ending here under SemiGlobal (free trailing
        // gap, the default) isolates the forward-propagation question: does
        // the still-untouched, never-reconverged 5..15 chain get its cached
        // nearest_fork updated once node 4 retroactively becomes a fork.
        let mut divergent = b("AAAAC");
        divergent.push(b'T'); // seed has 'G' here; this read has 'T'
        for _ in 0..3 {
            g.add_read(&divergent).unwrap();
        }
        assert_eq!(
            g.nodes[fork_node].out_edges.len(),
            2,
            "node {fork_node} should now be a genuine fork"
        );

        // Independent ground truth: walk backward from each descendant
        // through raw in_edges (exhaustive, no hop bound -- safe here since
        // the test graph is tiny), looking for the nearest node with 2+
        // out-edges. This is deliberately NOT the cache and NOT the
        // production interior-filter walk, so agreement is real evidence,
        // not the cache checking itself.
        fn ground_truth_nearest_fork(nodes: &[Node], start: usize) -> Option<(usize, usize)> {
            let mut cur = start;
            loop {
                if nodes[cur].in_edges.len() != 1 {
                    return None; // reached a reconvergence/source with no single predecessor
                }
                let pred = nodes[cur].in_edges[0];
                if nodes[pred].out_edges.len() >= 2 {
                    // `cur` is pred's direct child -- the arm entry point --
                    // not the node we started the walk from.
                    return Some((pred, cur));
                }
                cur = pred;
            }
        }

        for idx in 5..seed.len() {
            let expected = ground_truth_nearest_fork(&g.nodes, idx);
            assert_eq!(
                g.nodes[idx].nearest_fork, expected,
                "node {idx}: cached nearest_fork should match the independently-computed \
                 ground truth after node {fork_node} became a fork post-hoc"
            );
            assert_eq!(
                g.nodes[idx].nearest_fork,
                Some((fork_node, orig_child)),
                "node {idx}: expected the cache to have been updated by forward \
                 propagation to point at the newly-created fork"
            );
        }
    }

    /// The adversarial case that the first version of the test above
    /// actually hit (by accident, before its own construction was fixed):
    /// a single read both creates a brand-new fork AND, later in that exact
    /// same read's own traceback, reconverges back onto the pre-existing
    /// chain. This is exactly the timing hazard eager, inline propagation
    /// got wrong (a downstream node looked like an ordinary single-
    /// predecessor node at the moment the fork upstream was created, then
    /// gained a second in-edge moments later in the same read). Confirms
    /// the reconvergence node itself, and everything downstream of it, is
    /// correctly left untouched (still `None`, matching an independent
    /// ground-truth walk) rather than being stale-updated to point at the
    /// new fork.
    #[test]
    fn nearest_fork_same_read_reconvergence_is_not_stale_updated() {
        let seed = b("AAAACGTACGTACGTA");
        let cfg = PoaConfig {
            min_reads: 3,
            band_width: 0,
            adaptive_band: false,
            warn_on_long_unbanded: false,
            ..Default::default()
        };
        let mut g = PoaGraph::new(&seed, cfg).unwrap();
        for _ in 0..4 {
            g.add_read(&seed).unwrap();
        }
        let fork_node = 4;
        let orig_child = g.nodes[fork_node].out_edges[0].0;
        assert_eq!(orig_child, 5);

        // This divergent read substitutes one base at position 5, then
        // continues with the *same* suffix as the seed (seed[6..]) -- so
        // after the mismatch, its remaining bases genuinely re-Match nodes
        // 6..15 one-for-one, giving node 6 a real second in-edge from the
        // new mismatch node. Node 4 becomes a fork; node 6 becomes a true
        // reconvergence point; both effects come from this one read.
        let mut divergent = b("AAAAC");
        divergent.push(b'T'); // seed has 'G' here
        divergent.extend_from_slice(&seed[6..]);
        for _ in 0..3 {
            g.add_read(&divergent).unwrap();
        }
        assert_eq!(
            g.nodes[fork_node].out_edges.len(),
            2,
            "node {fork_node} should now be a genuine fork"
        );
        let reconv = 6usize;
        assert_eq!(
            g.nodes[reconv].in_edges.len(),
            2,
            "node {reconv} should have genuinely reconverged (2 predecessors) \
             within the same read that created the fork"
        );

        fn ground_truth_nearest_fork(nodes: &[Node], start: usize) -> Option<(usize, usize)> {
            let mut cur = start;
            loop {
                if nodes[cur].in_edges.len() != 1 {
                    return None;
                }
                let pred = nodes[cur].in_edges[0];
                if nodes[pred].out_edges.len() >= 2 {
                    return Some((pred, cur));
                }
                cur = pred;
            }
        }

        // Node 5 sits strictly between the fork and the reconvergence: it
        // should still get the new fork.
        assert_eq!(
            g.nodes[5].nearest_fork,
            Some((fork_node, orig_child)),
            "node 5 (between the fork and the reconvergence) should still be updated"
        );
        assert_eq!(
            g.nodes[5].nearest_fork,
            ground_truth_nearest_fork(&g.nodes, 5)
        );

        // Node 6 (the reconvergence point itself) and everything downstream
        // of it must NOT be stale-updated to the new fork -- they should
        // match ground truth, which is None past a real reconvergence.
        for idx in 6..seed.len() {
            let expected = ground_truth_nearest_fork(&g.nodes, idx);
            assert_eq!(
                expected, None,
                "sanity check on the ground-truth helper itself: node {idx} \
                 is past a real reconvergence, so ground truth must be None"
            );
            assert_eq!(
                g.nodes[idx].nearest_fork, expected,
                "node {idx}: must not be stale-updated across the reconvergence \
                 at node {reconv} -- this is the exact hazard deferred (end-of-read) \
                 propagation is meant to avoid"
            );
        }
    }
}

// ─── White-box tests: bypass edges (Phase 1 of bypass_edge_delete_rework) ────
//
// `PoaGraph.bypass_edges` is a private field written by `add_to_graph` but not
// yet read by any consumer (Phase 1 is purely additive; wiring it into
// `heaviest_path` is Phase 2). These same-module tests assert the field's own
// contents directly against hand-constructed read sets with known deletion
// patterns, since nothing observable in `consensus()` output depends on it
// yet. All use unbanded Global alignment over a deliberately non-repetitive
// seed so every read's Match/Delete op sequence -- and therefore the exact
// (from, to) node indices of the bypass edges -- is unambiguous and exactly
// predictable (a linear seed builds nodes 0..len in sequence order, so for a
// pure-deletion read every op references node index == seed position).
#[cfg(test)]
mod bypass_edge_tests {
    use super::*;

    fn b(s: &str) -> Vec<u8> {
        s.as_bytes().to_vec()
    }

    // 16 bp, non-repetitive: no adjacent duplicate bases near the deletion
    // sites used below, so each deletion has a single optimal placement.
    const SEED: &str = "ACGTGCATCGTAGCTA";

    fn cfg_global_unbanded() -> PoaConfig {
        PoaConfig {
            band_width: 0,
            adaptive_band: false,
            alignment_mode: AlignmentMode::Global,
            warn_on_long_unbanded: false,
            min_reads: 2,
            ..Default::default()
        }
    }

    /// Flatten `bypass_edges` into a sorted `(from, to, weight)` list for
    /// order-independent assertions.
    fn flatten(g: &PoaGraph) -> Vec<(usize, usize, i32)> {
        let mut out: Vec<(usize, usize, i32)> = g
            .bypass_edges
            .iter()
            .flat_map(|(&from, tos)| tos.iter().map(move |&(to, w)| (from, to, w)))
            .collect();
        out.sort_unstable();
        out
    }

    /// The `EdgeWeight` on the `from -> to` out-edge, or `None` if no such edge.
    fn out_edge(g: &PoaGraph, from: usize, to: usize) -> Option<EdgeWeight> {
        g.nodes[from]
            .out_edges
            .iter()
            .find(|&&(t, _)| t == to)
            .map(|&(_, ew)| ew)
    }

    /// Asserts the pure-bypass invariant across the whole graph: no `out_edge`
    /// anywhere carries any `deleted`-bucket weight (under pure bypass a
    /// deleting read never traverses an edge into a skipped node), and
    /// `edge_delete_reads` was never populated.
    fn assert_no_delete_bucket_or_delete_reads(g: &PoaGraph) {
        for (idx, nd) in g.nodes.iter().enumerate() {
            for &(to, ew) in &nd.out_edges {
                assert_eq!(
                    ew.deleted, 0,
                    "pure bypass: out-edge {idx}->{to} must carry no deleted-bucket weight, \
                     got {}",
                    ew.deleted
                );
            }
        }
        assert!(
            g.edge_delete_reads.is_empty(),
            "pure bypass: edge_delete_reads must never be populated"
        );
    }

    #[test]
    fn single_base_deletion_records_one_bypass_edge() {
        // Delete seed position 6 ('A', between 'C'@5 and 'T'@7, both distinct
        // from 'A' and from each other -> unambiguous single-base deletion).
        let mut g = PoaGraph::new(&b(SEED), cfg_global_unbanded()).unwrap();
        let read = b("ACGTGCTCGTAGCTA"); // SEED without index 6
        g.add_read(&read).unwrap();

        // Expected op sequence: Match(0..=5), Delete(6), Match(7..=15).
        // Entry predecessor of the run = node 5; clean-match resume onto the
        // existing node 7 -> one bypass edge 5->7, weight 1.
        assert_eq!(
            flatten(&g),
            vec![(5, 7, 1)],
            "one single-base deletion should record exactly one bypass edge 5->7 weight 1"
        );
        // (a) delete_count on the skipped node is still incremented -- this is
        // all mean_column_entropy/majority_frequency need.
        assert_eq!(
            g.nodes[6].delete_count, 1,
            "delete_count on the skipped node must still be incremented"
        );
        // Pure-bypass representation: the laundered through-edge is ABSENT.
        // Under the superseded dual-bookkeeping the deleting read advanced
        // through node 6 and its resume Match(7) inflated edge 6->7's matched
        // weight to 2 (seed + this read); under pure bypass 6->7 stays at the
        // seed's matched=1, and node 5's seed edge 5->6 gains no deleted-bucket
        // weight either.
        assert_eq!(
            out_edge(&g, 6, 7).map(|ew| ew.matched),
            Some(1),
            "edge 6->7 must NOT have gained matched weight from the deleting read \
             (that would be the laundering this rework removes)"
        );
        assert_no_delete_bucket_or_delete_reads(&g);
    }

    #[test]
    fn multi_base_deletion_run_records_one_spanning_bypass_edge() {
        // Delete seed positions 6,7 ('A','T'): "C A T C"@5..=8 -> "C C",
        // a clean 2-base deletion between the two 'C's.
        let mut g = PoaGraph::new(&b(SEED), cfg_global_unbanded()).unwrap();
        let read = b("ACGTGCCGTAGCTA"); // SEED without indices 6,7
        g.add_read(&read).unwrap();

        // Match(0..=5), Delete(6), Delete(7), Match(8..=15): a single bypass
        // edge spanning the whole run, 5 -> 8, not one edge per deleted node.
        assert_eq!(
            flatten(&g),
            vec![(5, 8, 1)],
            "a 2-base deletion run should record exactly one bypass edge 5->8 weight 1"
        );
        assert_eq!(g.nodes[6].delete_count, 1);
        assert_eq!(g.nodes[7].delete_count, 1);
        // Pure-bypass: resume edge 7->8 not inflated by the deleting read.
        assert_eq!(out_edge(&g, 7, 8).map(|ew| ew.matched), Some(1));
        assert_no_delete_bucket_or_delete_reads(&g);
    }

    #[test]
    fn deletion_run_spanning_a_preexisting_fork() {
        let mut g = PoaGraph::new(&b(SEED), cfg_global_unbanded()).unwrap();

        // First, an insertion read that gives node 7 ('T') a second out-edge
        // (to the first inserted node), making it a genuine fork. This read
        // has no Delete ops, so it records no bypass edge itself.
        let ins = b("ACGTGCATAACGTAGCTA"); // SEED with "AA" inserted after node 7
        g.add_read(&ins).unwrap();
        assert_eq!(
            g.nodes[7].out_edges.len(),
            2,
            "node 7 should now be a fork (matched arm to C@8, inserted arm to the new 'A')"
        );
        assert!(
            flatten(&g).is_empty(),
            "an insertion-only read must not record any bypass edge"
        );

        // Now a deletion read that deletes nodes 6,7,8 -- a run that spans the
        // fork at node 7. Added twice to also exercise weight incrementing.
        let del = b("ACGTGCGTAGCTA"); // SEED without indices 6,7,8
        g.add_read(&del).unwrap();
        g.add_read(&del).unwrap();

        // Match(0..=5), Delete(6,7,8), Match(9..=15): one bypass edge 5 -> 9,
        // weight 2 (two identical deletion reads). The fork at node 7 in the
        // middle of the deleted run does not fragment or misplace it.
        assert_eq!(
            flatten(&g),
            vec![(5, 9, 2)],
            "a deletion run spanning a fork should still record one bypass edge 5->9, \
             incremented to weight 2 across the two identical deletion reads"
        );
        assert_eq!(g.nodes[6].delete_count, 2);
        assert_eq!(g.nodes[7].delete_count, 2);
        assert_eq!(g.nodes[8].delete_count, 2);
        // Pure-bypass: resume edge 8->9 is matched=2 (the seed plus the
        // insertion read, both genuine matchers that traverse 8->9), NOT 4 --
        // the two DELETING reads did not launder their weight onto it. Under
        // the superseded dual-bookkeeping it would have been 4.
        assert_eq!(out_edge(&g, 8, 9).map(|ew| ew.matched), Some(2));
        assert_no_delete_bucket_or_delete_reads(&g);
    }

    #[test]
    fn leading_delete_run_records_no_bypass_edge() {
        // Leading delete run: read is a strict suffix of the seed. Under
        // Global alignment the seed's leading bases must be deleted to connect
        // to the source, so the read's ops begin with a Delete run -- one that
        // has no entry predecessor (`bypass_pending == Some(None)`). It
        // resolves (closes) at the first Match but records no bypass edge,
        // because there is no predecessor node to bypass *from*.
        let mut g = PoaGraph::new(&b(SEED), cfg_global_unbanded()).unwrap();
        let suffix = b("ATCGTAGCTA"); // SEED[6..16]
        g.add_read(&suffix).unwrap();
        assert!(
            flatten(&g).is_empty(),
            "a leading delete run (no predecessor to bypass from) must record no bypass edge"
        );
        // Confirm the leading Delete run genuinely occurred (so this really
        // exercises the Some(None) path, not a vacuous empty-ops case) while
        // the existing per-node delete bookkeeping still happens.
        assert_eq!(
            g.nodes[0].delete_count, 1,
            "leading deleted nodes are still counted by the existing bookkeeping"
        );
    }

    #[test]
    fn trailing_terminal_delete_run_records_no_bypass_edge() {
        // A trailing terminal Delete run (a Delete run with no following
        // Match/Insert to resume at) does not arise from `align()` over a
        // linear seed -- its terminal-cell search gives a free trailing gap,
        // so the optimal traceback simply ends at the last matched node rather
        // than deleting onward. To exercise the end-of-ops drop path directly
        // and deterministically, call `add_to_graph` with a hand-built op
        // sequence that ends in a Delete run (a shape a real traceback can
        // produce in a branching graph): Match(0), Match(1), Delete(2),
        // Delete(3). The pending bypass run must be dropped, never recorded as
        // an edge to a nonexistent resume node.
        let mut g = PoaGraph::new(&b(SEED), cfg_global_unbanded()).unwrap();
        let (_topo, rank_of) = topological_order(&g.nodes);
        let ops = vec![
            AlignOp::Match(0),
            AlignOp::Match(1),
            AlignOp::Delete(2),
            AlignOp::Delete(3),
        ];
        let query = b("AC"); // matches nodes 0 ('A') and 1 ('C')
        add_to_graph(
            &mut g.nodes,
            &mut g.edge_reads,
            &mut g.bypass_edges,
            &mut g.fork_arm_index,
            &rank_of,
            &query,
            &ops,
            1,
        );
        assert!(
            flatten(&g).is_empty(),
            "a trailing terminal delete run must not create a bypass edge to nowhere"
        );
        // The per-node delete bookkeeping still runs for the terminal deletes
        // -- only the (correctly absent) bypass edge differs.
        assert_eq!(g.nodes[2].delete_count, 1);
        assert_eq!(g.nodes[3].delete_count, 1);
        // Pure-bypass: the terminal Delete run created no edges into/out of the
        // skipped nodes.
        assert_no_delete_bucket_or_delete_reads(&g);
    }

    #[test]
    fn no_deletions_records_no_bypass_edges() {
        // A pure-match read set never touches bypass_edges at all.
        let mut g = PoaGraph::new(&b(SEED), cfg_global_unbanded()).unwrap();
        g.add_read(&b(SEED)).unwrap();
        g.add_read(&b(SEED)).unwrap();
        assert!(
            flatten(&g).is_empty(),
            "reads with no Delete ops must leave bypass_edges empty"
        );
    }
}

// ─── White-box tests: validate_and_merge_groups structure-corroborated split ─
//
// Phase 3 of design/bypass_edge_delete_rework.md added a "clean distinguishing
// bubble" requirement to `validate_and_merge_groups` (a private fn), so a
// single allele whose deletion-heavy reads `phasing_groups` split into
// narrowly length-separated sub-groups is merged back, while genuine alleles
// with a real arm-choice difference stay split. These tests exercise the
// decision function directly with hand-built `groups`/`bubbles`/`edge_reads`
// -- deterministic, and guarding BOTH directions of the discriminator (the
// full pbsim-driven integration case is the `multi_skew_cag20_40` bench
// scenario, which black-box reproduction of the exact over-split is too
// error-model-dependent to capture reliably in a unit test).
#[cfg(test)]
mod validate_and_merge_tests {
    use super::*;

    fn read_of_len(n: usize) -> Vec<u8> {
        vec![b'A'; n]
    }

    /// Two length-separated sub-groups that take the SAME majority arm at every
    /// bubble (structurally indistinguishable -- the multi_skew_cag20_40 shape:
    /// one allele split by deletion noise) must merge back into one allele.
    /// Without the fix, the length-separated sub-group is confirmed as its own
    /// allele and this returns 3 groups.
    #[test]
    fn merges_structurally_indistinct_length_split() {
        // reads 0..10 = short allele (~100bp); 10..15 = long sub-group A
        // (~316bp); 15..19 = long sub-group B (~300bp, narrowly length-
        // separated from A -- a 16bp median gap, like the real over-split).
        let mut reads = Vec::new();
        for _ in 0..10 {
            reads.push(read_of_len(100));
        }
        for k in 0..5 {
            reads.push(read_of_len(314 + k)); // 314..=318, median 316
        }
        for k in 0..4 {
            reads.push(read_of_len(298 + k)); // 298..=301, median 300
        }
        let g0: Vec<usize> = (0..10).collect();
        let sub_a: Vec<usize> = (10..15).collect();
        let sub_b: Vec<usize> = (15..19).collect();

        // Two structural bubbles (require_length_separation gate). Synthetic
        // node ids; arm 0 = "short" arm, arm 1 = "long" arm.
        let bubbles = vec![
            (0usize, vec![1usize, 2usize]),
            (3usize, vec![4usize, 5usize]),
        ];
        let mut edge_reads: HashMap<(usize, usize), Vec<u32>> = HashMap::new();
        edge_reads.insert((0, 1), g0.iter().map(|&r| r as u32).collect());
        edge_reads.insert((3, 4), g0.iter().map(|&r| r as u32).collect());
        // Both long sub-groups take arm 1 at both bubbles -> no distinguishing
        // bubble between them.
        let long: Vec<u32> = sub_a.iter().chain(&sub_b).map(|&r| r as u32).collect();
        edge_reads.insert((0, 2), long.clone());
        edge_reads.insert((3, 5), long);

        let out = validate_and_merge_groups(
            vec![g0.clone(), sub_a.clone(), sub_b.clone()],
            &reads,
            3,
            &bubbles,
            &edge_reads,
        );
        assert_eq!(
            out.len(),
            2,
            "structurally-indistinct length-split sub-groups must merge to 2 alleles, got {}",
            out.len()
        );
        // Every long-allele read (both sub-groups) ends up together.
        let long_group = out
            .iter()
            .find(|g| g.contains(&10))
            .expect("a group containing the long allele");
        for r in sub_a.iter().chain(&sub_b) {
            assert!(
                long_group.contains(r),
                "long-allele read {r} must be in the single merged long group"
            );
        }
    }

    /// Companion (the over-merge guard): two genuine alleles that DO have a
    /// clean distinguishing bubble (different majority arms) must stay split --
    /// the fix narrows splitting, so it must not collapse real alleles.
    #[test]
    fn keeps_structurally_distinct_length_split() {
        let mut reads = Vec::new();
        for _ in 0..10 {
            reads.push(read_of_len(100)); // allele X, arm 0
        }
        for _ in 0..8 {
            reads.push(read_of_len(200)); // allele Y, arm 1
        }
        let gx: Vec<usize> = (0..10).collect();
        let gy: Vec<usize> = (10..18).collect();
        let bubbles = vec![
            (0usize, vec![1usize, 2usize]),
            (3usize, vec![4usize, 5usize]),
        ];
        let mut edge_reads: HashMap<(usize, usize), Vec<u32>> = HashMap::new();
        edge_reads.insert((0, 1), gx.iter().map(|&r| r as u32).collect());
        edge_reads.insert((3, 4), gx.iter().map(|&r| r as u32).collect());
        edge_reads.insert((0, 2), gy.iter().map(|&r| r as u32).collect());
        edge_reads.insert((3, 5), gy.iter().map(|&r| r as u32).collect());

        let out = validate_and_merge_groups(vec![gx, gy], &reads, 3, &bubbles, &edge_reads);
        assert_eq!(
            out.len(),
            2,
            "two genuine alleles with a clean distinguishing bubble must stay split, got {}",
            out.len()
        );
    }
}

// ─── White-box tests: BandTooNarrow + align_with_retry ─────────────────────
//
// `align` and `align_with_retry` are private, so directly confirming the
// exact threshold at which a too-narrow band now returns `BandTooNarrow`
// (rather than the old silent empty-ops collapse), and that the retry
// wrapper recovers a correct alignment, requires a same-module test.
//
// Fixture: a 100bp seed and a 220bp query built by appending 120bp of new
// content to the seed's own sequence -- deliberately non-repetitive (unlike
// `tests/adaptive_band_collapse.rs`'s real GAA-repeat fixture) so the
// expected alignment shape (100 Match + 120 Insert, 0 Delete) is
// unambiguous and exactly predictable, isolating the band/retry mechanism
// itself from any repeat-driven alignment ambiguity.
#[cfg(test)]
mod band_too_narrow_tests {
    use super::*;

    fn fixture() -> (Vec<u8>, Vec<u8>) {
        let seed: Vec<u8> = "ACGTACGTCG".repeat(10).into_bytes(); // 100bp
        let mut query = seed.clone();
        query.extend_from_slice(&"TGCA".repeat(30).into_bytes()); // +120bp
        (seed, query)
    }

    fn align_direct(
        g: &PoaGraph,
        topo: &[usize],
        rank_of: &[usize],
        spine: &[(usize, u8, i32)],
        query: &[u8],
        band_width: usize,
        scratch: &mut AlignScratch,
    ) -> Result<Vec<AlignOp>, PoaError> {
        let cfg = PoaConfig {
            band_width,
            adaptive_band: false,
            alignment_mode: AlignmentMode::SemiGlobal,
            ..PoaConfig::default()
        };
        align(&g.nodes, topo, rank_of, spine, query, &cfg, scratch, &[])
    }

    /// Directly confirms the exact band-width threshold: below it, `align`
    /// now returns `Err(BandTooNarrow)` (never an empty, silently-wrong
    /// `Ok(vec![])` -- the bug `tests/adaptive_band_collapse.rs` guards
    /// against); at or above it, `align` succeeds directly with no retry
    /// needed. Values confirmed empirically when this test was written, not
    /// assumed: 100 (exactly the seed's own length, plausible-looking but
    /// still 120bp short of the true diagonal shift) still errors; 119
    /// (one short of the full 120bp gap) already succeeds, because the
    /// windowing here centres with a small amount of slack, not because
    /// the boundary is off by one in the traceback itself.
    #[test]
    fn band_too_narrow_returned_below_threshold_ok_at_and_above() {
        let (seed, query) = fixture();
        let g = PoaGraph::new(&seed, PoaConfig::default()).unwrap();
        let (topo, rank_of) = topological_order(&g.nodes);
        let spine = heaviest_path(&g.nodes, &topo, &rank_of, &g.bypass_edges);
        let mut scratch = AlignScratch::new();

        for &band_width in &[10usize, 50, 100] {
            let result = align_direct(
                &g,
                &topo,
                &rank_of,
                &spine,
                &query,
                band_width,
                &mut scratch,
            );
            match result {
                Err(PoaError::BandTooNarrow {
                    configured,
                    required,
                }) => {
                    assert_eq!(configured, band_width);
                    assert!(
                        required > band_width,
                        "required ({required}) should exceed the too-narrow configured \
                         width ({band_width})"
                    );
                }
                other => panic!(
                    "band_width={band_width} expected BandTooNarrow, got {:?}",
                    other.map(|ops| ops.len())
                ),
            }
        }

        for &band_width in &[119usize, 120, 121, 150] {
            let result = align_direct(
                &g,
                &topo,
                &rank_of,
                &spine,
                &query,
                band_width,
                &mut scratch,
            );
            match result {
                Ok(ops) => assert_eq!(
                    ops.len(),
                    query.len(),
                    "band_width={band_width}: expected one op per query base \
                     (100 Match + 120 Insert), got {} ops",
                    ops.len()
                ),
                Err(e) => panic!("band_width={band_width} expected Ok, got {e:?}"),
            }
        }
    }

    /// Confirms `align_with_retry` recovers a fully correct alignment from
    /// every narrow starting width tested above, transparently -- the
    /// actual fix for callers, not just the diagnostic in the test above.
    #[test]
    fn align_with_retry_recovers_correct_alignment_from_any_narrow_start() {
        let (seed, query) = fixture();
        let g = PoaGraph::new(&seed, PoaConfig::default()).unwrap();
        let (topo, rank_of) = topological_order(&g.nodes);
        let spine = heaviest_path(&g.nodes, &topo, &rank_of, &g.bypass_edges);
        let mut scratch = AlignScratch::new();

        for &band_width in &[10usize, 50, 100] {
            let cfg = PoaConfig {
                band_width,
                adaptive_band: false,
                alignment_mode: AlignmentMode::SemiGlobal,
                ..PoaConfig::default()
            };
            let (ops, retried) = align_with_retry(
                &g.nodes,
                &topo,
                &rank_of,
                &spine,
                &query,
                &cfg,
                &mut scratch,
                &[],
            )
            .unwrap_or_else(|e| panic!("band_width={band_width}: retry should recover, got {e:?}"));
            assert!(
                retried,
                "band_width={band_width}: expected align_with_retry to report that pass 1 \
                 needed widening"
            );

            let n_match = ops
                .iter()
                .filter(|o| matches!(o, AlignOp::Match(_)))
                .count();
            let n_insert = ops
                .iter()
                .filter(|o| matches!(o, AlignOp::Insert(_)))
                .count();
            let n_delete = ops
                .iter()
                .filter(|o| matches!(o, AlignOp::Delete(_)))
                .count();
            assert_eq!(
                (n_match, n_insert, n_delete),
                (100, 120, 0),
                "band_width={band_width}: expected the retry to recover the exact \
                 alignment shape (100 Match + 120 Insert, 0 Delete), got \
                 Match={n_match} Insert={n_insert} Delete={n_delete}"
            );
        }
    }
}
