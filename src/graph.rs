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
#[inline]
fn edge_weight(nodes: &[Node], from: usize, to: usize) -> i32 {
    nodes[from]
        .out_edges
        .iter()
        .find(|&&(t, _)| t == to)
        .map(|&(_, w)| w)
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
    /// Reusable scratch for lock-window tables; avoids per-call heap allocation.
    align_scratch: AlignScratch,
    /// Minimizer index over the cached spine sequence; rebuilt whenever
    /// `cached_spine` is refreshed.  Maps k-mer hash → spine rank (index into
    /// `cached_spine`).  Only hashes that appear exactly once in the spine are
    /// stored — duplicate k-mers cannot serve as unambiguous anchors.
    spine_mers: HashMap<u64, u32>,
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

/// Half-width of the DP band applied to spine (non-bubble) nodes.
/// Minimum SPINE_MARGIN when nothing else constrains it (no band specified).
/// Covers ±SPINE_MARGIN_MIN query positions for unbanded or wide-band configs.
const SPINE_MARGIN_MIN: usize = 50;

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
    // Uses the same adaptive formula as main (b + f*L), floored at SPINE_MARGIN_MIN
    // only when band_width=0 and adaptive_band=false (unbanded / no guidance).
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
        SPINE_MARGIN_MIN
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
                    if bj < l && node_base == query[bj] {
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
        None => (0, State::M),
    };

    // ── Traceback ─────────────────────────────────────────────────────────────
    let mut ops: Vec<AlignOp> = Vec::with_capacity(l + n / 4);
    let mut t = best_t;
    let mut j = l;
    let mut cur_state = best_state;

    loop {
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
    let mut longest_bubble_span = 0usize;
    for nd in nodes {
        let qualifying: Vec<(usize, i32)> = nd
            .out_edges
            .iter()
            .filter(|&&(_, w)| w >= threshold)
            .map(|&(to, w)| (to, w))
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
            let arms: Vec<usize> = nodes[entry_node]
                .out_edges
                .iter()
                .filter(|&&(_, w)| w >= threshold)
                .map(|&(to, _)| to)
                .collect();

            if arms.len() < 2 {
                return None;
            }

            // An arm whose start node already has multiple in-edges is a direct edge to
            // the bubble exit (arm span 0). Otherwise measure the single-successor chain.
            let max_span = arms
                .iter()
                .map(|&start| {
                    if nodes[start].in_edges.len() > 1 {
                        0
                    } else {
                        materialize_arm_len(nodes, start, max_check)
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
                let len = if nodes[a].in_edges.len() > 1 {
                    0
                } else {
                    materialize_arm_len(nodes, a, phasing_bubble_min_span)
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

/// Groups reads by arm-choice compatibility across all structural bubbles.
///
/// For each structural bubble, edge_reads tells us which reads entered each arm.
/// Two reads are in the same haplotype group when they agree on every bubble
/// where both have a recorded arm choice. Reads with no arm assignment at any
/// structural bubble (pre-dating all bubbles, or not spanning them) are folded
/// into the largest group. Groups below min_reads are also folded into the largest.
fn phasing_groups(
    edge_reads: &HashMap<(usize, usize), Vec<u32>>,
    bubbles: &[(usize, Vec<usize>)],
    n_reads: usize,
    min_reads: usize,
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

    // Union-find over assigned reads: merge reads whose arm choices never conflict.
    let n = assigned.len();
    let mut parent: Vec<usize> = (0..n).collect();

    for i in 0..n {
        for j in (i + 1)..n {
            let ri = assigned[i];
            let rj = assigned[j];
            let compatible = (0..n_bubbles).all(|b| match (sig[ri][b], sig[rj][b]) {
                (Some(a), Some(bv)) => a == bv,
                _ => true,
            });
            if compatible {
                let mut pi = i;
                while parent[pi] != pi {
                    pi = parent[pi];
                }
                let mut pj = j;
                while parent[pj] != pj {
                    pj = parent[pj];
                }
                if pi != pj {
                    parent[pj] = pi;
                }
            }
        }
    }

    // Path compression.
    for i in 0..n {
        let mut root = i;
        while parent[root] != root {
            root = parent[root];
        }
        parent[i] = root;
    }

    // Collect groups.
    let mut group_map: HashMap<usize, Vec<usize>> = HashMap::new();
    for (slot, &r) in assigned.iter().enumerate() {
        group_map.entry(parent[slot]).or_default().push(r);
    }

    let mut groups: Vec<Vec<usize>> = group_map.into_values().collect();
    groups.sort_unstable_by_key(|g| std::cmp::Reverse(g.len()));

    if groups.is_empty() {
        groups.push(Vec::new());
    }

    // Fold groups below min_reads and unassigned reads into the largest group.
    let mut i = 1;
    while i < groups.len() {
        if groups[i].len() < min_reads {
            let g = groups.remove(i);
            groups[0].extend(g);
        } else {
            i += 1;
        }
    }
    groups[0].extend(unassigned);
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
            align_scratch: AlignScratch::new(),
            spine_mers: HashMap::new(),
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
        // consensus always recomputes the heaviest path from scratch.  A stale
        // spine affects both alignment speed and anchor quality: shifted k-mers
        // built from an outdated spine can produce wrong anchor chains.
        let reads_since_update = self.n_reads.saturating_sub(self.spine_updated_at);
        if self.cached_spine.is_empty() || reads_since_update >= self.spine_interval {
            let new_spine = heaviest_path(&self.nodes, &topo, &rank_of);
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

        let ops = align(
            &self.nodes,
            &topo,
            &rank_of,
            &self.cached_spine,
            read,
            &self.config,
            &mut self.align_scratch,
            &anchors,
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
        let min_cov = self.min_coverage();

        let filtered: Vec<(usize, u8, i32)> = match self.config.consensus_mode {
            ConsensusMode::HeaviestPath => {
                let path = heaviest_path(&self.nodes, &topo, &rank_of);

                // Step 1: per-node support. A node is adequately supported
                // if its Match coverage clears the global threshold, OR it
                // is the locally dominant arm at its bubble: heaviest_path
                // already picked it as the winning predecessor edge, and its
                // own Match coverage clears a majority of its *own* bubble's
                // vote total, even though global coverage falls short
                // because other reads diverged at an earlier bubble (so
                // fewer than min_cov reads even reach this point). Gated on
                // the predecessor having 2+ out-edges (a real fork) so a
                // plain low-coverage unbranched run -- a true minority
                // extension, not a bubble arm -- is still governed by the
                // global check alone. Coverage (not edge weight) is compared
                // against the local threshold: edge weight includes reads
                // that Delete through the node (skip it without adopting its
                // base), which is evidence *against* keeping the base.
                let supported: Vec<bool> = path
                    .iter()
                    .enumerate()
                    .map(|(i, &(node_idx, _, _))| {
                        if self.nodes[node_idx].coverage >= min_cov {
                            return true;
                        }
                        let Some(&(pred_idx, _, _)) = i.checked_sub(1).and_then(|j| path.get(j))
                        else {
                            return false;
                        };
                        if self.nodes[pred_idx].out_edges.len() < 2 {
                            return false;
                        }
                        let local_total: i32 = self.nodes[pred_idx]
                            .out_edges
                            .iter()
                            .map(|&(_, ew)| ew)
                            .sum();
                        let local_min_cov = coverage_threshold(
                            local_total.max(0) as usize,
                            self.config.min_coverage_fraction,
                        );
                        self.nodes[node_idx].coverage as i32 >= local_min_cov as i32
                    })
                    .collect();

                let start = supported.iter().position(|&s| s).unwrap_or(0);
                let end = supported
                    .iter()
                    .rposition(|&s| s)
                    .map(|i| i + 1)
                    .unwrap_or(path.len());

                let (effective_path, effective_supported) = if start < end {
                    (&path[start..end], &supported[start..end])
                } else {
                    (&path[..], &supported[..])
                };

                effective_path
                    .iter()
                    .zip(effective_supported.iter())
                    .filter(|&(_, &s)| s)
                    .map(|(&t, _)| t)
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
        let ops = align(
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
        let spine = heaviest_path(&self.nodes, &topo, &rank_of);
        let ops = align(
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

    /// Return a snapshot of the graph topology for visualization and inspection.
    ///
    /// Nodes are in topological order; `edges` use topological ranks as source
    /// and target identifiers.  `spine_ranks` lists the ranks of nodes on the
    /// heaviest-path (consensus) spine.
    pub fn graph_topology(&self) -> crate::types::GraphTopology {
        use crate::types::{GraphEdgeInfo, GraphNodeInfo, GraphTopology};

        let (topo, rank_of) = topological_order(&self.nodes);
        let spine = heaviest_path(&self.nodes, &topo, &rank_of);

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

        let edges: Vec<GraphEdgeInfo> = topo
            .iter()
            .enumerate()
            .flat_map(|(from_rank, &node_idx)| {
                self.nodes[node_idx]
                    .out_edges
                    .iter()
                    .map(move |&(succ_idx, weight)| (from_rank, succ_idx, weight))
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
    ///
    /// Structural bubbles (arm span ≥ phasing_bubble_min_span) are phased first using
    /// cross-bubble compatibility grouping. If no structural bubbles are found and the
    /// alignment was banded, the graph is rebuilt with unbanded alignment so that large
    /// allele-length differences create a visible Insert-arm bubble; structural phasing
    /// then retries on the unbanded graph. If still no structural bubble, falls back to
    /// single-best SNP bubble partitioning for substitution haplotypes.
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
            phasing_groups(
                &self.edge_reads,
                &structural,
                self.n_reads,
                self.config.min_reads,
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
