# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Fixed

- **Follow-up to the multi-allele contamination fix below: the residual
  15.4% (8/52) misassignment rate reported after that fix was itself two
  more targeted, fixable bugs in `phasing_groups`/`validate_and_merge_groups`,
  not the periodic-alignment-ambiguity class (Known Bugs #3/#4/#9) it was
  provisionally attributed to.** Traced by reproducing the exact ground-truth
  ID tracing methodology on `multi_gaa30_100` and identifying each
  misassigned read's own recorded arm choice at every structural bubble.
  First finding, before any new fix: **the original 15.4%/8-read figure
  itself was wrong**, an off-by-construction bug in the *investigation*
  tooling, not the library -- `Consensus::read_indices` are indices into
  `PoaGraph`'s internal read order (seed = index 0, then `add_read` calls in
  external order skipping the seed), not the external `reads` slice passed
  to `consensus_multi`; the probe script used for the original trace treated
  them as external indices directly. Recomputing with the correct internal
  ↔ external mapping gives corrected baseline numbers: **32.7% (17/52)
  before the earlier fix, 11.5% (6/52) after it** -- still a large real
  improvement, just not the exact figures previously reported.
  - **(A) `phasing_groups` dumped every read it couldn't uniquely cluster
    (a "bridge candidate," compatible with 2+ clusters at once) into
    whichever cluster was globally largest, even when that read was
    provably *incompatible* with that cluster's own arm-choice template.**
    A bridge candidate's ambiguity means its bubble evidence didn't pick a
    *unique* cluster, not that the read is uninformative -- confirmed on
    `multi_gaa30_100`: 5 of the then-6 misassigned reads were bridge
    candidates whose compatible-cluster set never included the cluster they
    got dumped into. Fixed by resolving each bridge candidate (and each
    fully-`None`/uninformative read) against read length instead of cluster
    size: length is a signal largely orthogonal to bubble topology (a
    structural bubble exists *because* of a length difference), so it
    reliably breaks a tie the bubble evidence itself couldn't. Gated on a
    plausibility check (`PLAUSIBLE_LEN_FRACTION = 0.85` of the shortest
    well-supported cluster's median) before trusting a read's own length at
    all -- an unconditional first version of this fix regressed immediately:
    partial (non-spanning) reads have a truncated length reflecting how much
    of the flank+repeat they happened to cover, not their true allele's
    length, and several such reads clustered together purely by coincidental
    shortness, fabricating a spurious 3rd allele out of reads that were all
    genuinely the same true allele. A read that fails the plausibility check
    falls back to its largest *bubble-compatible* cluster by member count
    (bridge candidates) or the single largest cluster overall (fully
    uninformative reads) -- i.e. close to the previous behavior, but only
    for reads where a length comparison would itself be unreliable.
  - **(B) `validate_and_merge_groups` merged a whole rejected candidate group
    into its nearest confirmed group by the candidate's *aggregate* median
    length, which could still drag a genuinely different-origin read along
    with the majority of its own (accidentally mixed) cluster.** Confirmed
    on the one read (of 6) that fix (A) did not touch: 3 reads (2 true-short,
    1 true-long) shared byte-identical arm-choice signatures across all 3
    structural bubbles by coincidence -- a genuine per-read alignment tie,
    the one case in this investigation that really does match the
    periodic-alignment-ambiguity family -- so `phasing_groups` correctly (by
    its own logic) clustered all 3 together as one indivisible group. That
    group's aggregate median length was dominated by its 2-read majority and
    read as "short," so the whole group -- including the 1 true-long read --
    merged into the short confirmed group. Fixed by routing each member of a
    rejected candidate individually to whichever confirmed group its own
    length is nearest to, rather than moving the candidate as one block; safe
    specifically because the confirmed groups being compared against are, by
    construction, already length-separated from each other (that is what
    qualified them as confirmed), so per-read length comparison against them
    is not diluted by the same aggregation problem as the rejected
    candidate's own median.
  - **Net result on `multi_gaa30_100`:** misassignment rate goes from a
    corrected 11.5% (6/52, after the earlier fix alone) to **0.0% (0/52)**
    after both (A) and (B) -- every read in the scenario, reproducibly.
    Fix (A) alone accounted for 3.8% (2/52) as an intermediate checkpoint
    (with a regression to a spurious 3rd allele from an unconditional first
    attempt, caught and corrected before finalizing -- see above); (A) with
    its plausibility guard alone reached 1.9% (1/52); (B) closed the last
    read. `bench/validate.py` improves from 18/20 to **19/20**:
    `multi_gaa30_100` itself was already passing after the earlier fix and
    remains passing (Δ-2/Δ+1, unchanged); `multi_skew_cag20_40` -- previously
    reported as a FAIL attributed to "unrelated per-partition consensus
    noise" -- flips to a clean OK (Δ+0/Δ+0). That earlier attribution was
    itself reassessed: re-tracing `multi_skew_cag20_40`'s read provenance
    after this round's fix still shows a real, nonzero 12.8% (5/39)
    misassignment rate, smaller than `multi_gaa30_100`'s pre-fix rate but not
    zero -- the *consensus sequence* recovered correctly regardless, because
    the contaminating reads stayed a minority within their misassigned
    group's own heaviest-path vote. This is reported plainly rather than
    claimed as a second full fix: bench-level pass and zero read-level
    misassignment are related but not equivalent, and only the former is
    confirmed clean here. The remaining single-allele `cag50_d20` FAIL is
    unrelated (not a multi-allele scenario) and unaffected either way.
    `bench/compare_callers.py` stays 16/16 (no multi-allele scenarios in that
    suite). Full `cargo test --release --features cli` (192 lib tests, all
    integration binaries, doctests) green, no regressions, including the two
    Phase 2 tests (`partition_reads_by_bubble_excludes_delete_only_reads`,
    `bubble_site_arm_read_counts_excludes_delete_only_reads`) and both new
    tests from the earlier fix (`structural_phasing_no_contamination_on_noisy_periodic_repeat`,
    `structural_bubble_phasing_subtle_length_delta_not_rejected`). `cargo
    clippy --all-features --tests` and `cargo fmt --check` clean against
    this round's changes. `phasing_groups`'s signature gained a `reads: &[Vec<u8>]`
    parameter (it is a private function; no public API change).
    `PoaConfig`'s defaults are unchanged.

- **Multi-allele read-partition contamination on periodic (homogeneous) tandem
  repeats: `consensus_multi` could silently misassign a large fraction of reads
  to the wrong allele group.** Root-caused via ground-truth read tracing on real
  GAA×30/GAA×100 HiFi data (`multi_gaa30_100`): ~29% of reads (15 of 52) ended
  up in the wrong output allele's group. Two independent, layered fixes, applied
  in a mandated order so each one's own contribution could be measured before
  stacking the next:
  - **(a) `find_structural_bubbles`'s arm-span measurement is now noise-tolerant.**
    `materialize_arm_len` walked a strictly unbranched chain and stopped at the
    *first* fork or reconvergence, no matter how small -- a single read's indel
    error. In a periodic repeat, ordinary per-read indel noise (the dominant
    ONT/HiFi error mode in these regions -- matching-vs-inserting at any one
    column is close to a coin flip) scatters many such tiny disruptions
    throughout the entire length-differentiating region, so no arm ever
    measured long enough to clear `phasing_bubble_min_span` and the whole
    locus fell through to the single-column SNP-bubble fallback, which cannot
    reliably track a read's true (much longer) haplotype identity. Added
    `materialize_arm_len_tolerant` plus `real_in_edge_count`: a node's
    in-edges (or a fork's out-edges) only count as a "real" reconvergence/fork
    if the edge's own weight clears the same significance threshold already
    used elsewhere (`ceil(n_reads * min_allele_freq)`, floored at 1) --
    concretely, "minor branching" means *one or more* out-edges/in-edges whose
    weight is below that population-scaled floor, which by construction can
    never itself be a second true haplotype (a true second allele's arm must,
    by definition, already have cleared `min_allele_freq` to be worth phasing
    at all). The walk passes straight through any such minor branching and
    keeps measuring the arm; it still stops the instant *two or more*
    sufficiently-weighted arms/in-edges appear, so a genuine second
    structural bubble immediately downstream is never swallowed into the
    first arm's span. `collect_bubble_sites`'s public `BubbleSite.is_structural`
    field now uses the same tolerant measurement for consistency with
    `find_structural_bubbles`'s own classification -- **this is a computed-value
    behavior change on public API surface**: a site that was `is_structural:
    false` purely because of scattered per-read noise fragmenting its span may
    now report `true`. No other `materialize_arm_len` call site
    (`compute_stats`'s `longest_bubble_span`, the public `bubble_arm_lengths()`)
    was touched -- both remain general diagnostics, not consensus-affecting
    decisions, and are out of scope here.
  - **`phasing_groups`'s pairwise union-find had a latent wildcard-bridging bug,
    surfaced (not introduced) by (a) finding more structural bubbles than
    before.** A read with no recorded arm choice at some bubble (a wildcard,
    `None`) union-found its way into merging two clusters that a fully-informed
    read would never have judged compatible, because plain union-find only
    checks pairwise compatibility, not compatibility against a cluster's *entire*
    accumulated signature. Rewrote as informativeness-ordered incremental
    clustering: reads are processed most-fully-specified first; a read only
    joins an existing cluster if it is compatible with *every* bubble that
    cluster has resolved so far, joins nothing (starts a new cluster) if
    compatible with none, and is deferred to `bridge_candidates` (folded into
    the largest final group) if compatible with more than one -- never silently
    picking one arbitrarily. Verified this still correctly separates a nested
    3-allele case (`structural_bubble_phasing_three_alleles`, unchanged) while
    no longer bridging through partial-information reads.
  - **(b) `validate_and_merge_groups`: a read-length bimodality safety net for
    whatever the structural-bubble path still hands off.** After (a) also fixed
    the underlying union-find bug, provisional groups from `phasing_groups`
    could still include spurious small clusters (noise, not real minority
    alleles) that individually cleared the read-count floor. Before trusting a
    provisional split, each candidate group's read lengths are compared against
    already-confirmed groups' lengths via a median-gap-vs-spread test
    (`length_separated`: median gap >= 8bp *and* >= 3x the pooled median
    absolute deviation, MAD floored at 3bp so near-zero-noise clean read sets
    don't produce a degenerate zero threshold). A rejected candidate is folded
    into whichever confirmed group its own median length is closest to, rather
    than discarded. Two points tuned empirically, not assumed correct by design,
    documented in code:
    - The length-separation check only runs when the locus has **2 or more**
      structural bubbles, not universally. A single structural bubble can
      legitimately represent a same-length substitution-type difference (e.g.
      `locked_arm_deep_bubble_alleles_lost`'s G×30-vs-A×30, both 50bp) with no
      selection risk to guard against; the coincidental-tie risk this check
      exists for is specifically the multiple-comparisons exposure that appears
      once 2+ independent bubble sites are involved.
    - The significance bar for a candidate group is the plain `min_reads` floor,
      not a `min_allele_freq`-derived population-wide threshold. The stricter
      threshold double-penalizes genuine minority alleles, because
      `phasing_groups`'s own cross-bubble-consistency requirement already
      attrites a minority cluster's size below its raw single-bubble arm
      weight; re-applying the un-attrited population-wide bar on top of that
      wrongly rejected a confirmed 5-read CAG×40 minority cluster
      (`multi_skew_cag20_40`) that should have passed.
    - **This check is deliberately scoped to the structural-bubble/`phasing_groups`
      path only, not the same-length SNP-bubble fallback
      (`partition_reads_by_bubble`)** -- genuine same-length substitution
      haplotypes (`structural_bubble_phasing_ignores_snp_bubbles`,
      `CATCATCAT` vs `CATCGTCAT`) legitimately have zero length difference by
      construction, so a length check on that path would reject every correct
      call.
  - **Validation:** ground-truth read-provenance tracing on `multi_gaa30_100`
    (extracted reads' true allele label, from the simulator, cross-referenced
    against `Consensus.read_indices`) shows the misassignment rate drop from
    **28.8% (15/52) before this fix to 15.4% (8/52) after** -- a substantial
    improvement, not a complete elimination; the residual is consistent with
    Known Bug #3/#4's separately-tracked periodic-alignment ambiguity, not this
    contamination mechanism. `bench/validate.py`: `multi_gaa30_100` itself moves
    from FAIL to a clean OK (100->98 units, Δ-2; 30->31 units, Δ+1, both within
    tolerance), for an overall 18/20 (unchanged scenario count from the
    preceding harness fix, since that fix already resolved `sv_cag20_out60`
    separately -- see that entry below). `multi_skew_cag20_40` was investigated
    as a residual FAIL (20->22, Δ+2 on the majority allele) and confirmed, via
    the same ground-truth tracing method, to have **zero cross-contamination**
    post-fix -- the residual is ordinary per-partition consensus noise (the
    separately-tracked seed-sensitivity class of issue), not this bug, and is
    reported here rather than swept aside. `bench/compare_callers.py`: 16/16,
    unchanged (none of its scenarios exercise multi-allele partitioning).
    New regression tests: `tests/multi_allele_periodic_phasing.rs`
    (`structural_phasing_no_contamination_on_noisy_periodic_repeat`, a
    synthetic GAA×15-vs-GAA×40 read set with deterministic scattered indel
    noise -- confirmed to fail against the pre-fix code with ~30% per-group
    cross-contamination, matching the real-data magnitude, and pass after) and
    `structural_bubble_phasing_subtle_length_delta_not_rejected` in
    `src/tests.rs` (a genuine but subtle, 12bp, two-bubble length delta,
    confirming (b)'s bimodality check does not produce a false negative and
    silently downgrade a real heterozygous call to single-allele output). Full
    `cargo test --release --features cli` (192 lib tests, all integration
    binaries, doctests) green, including the two multi-allele partitioning
    tests from Phase 2 below (`partition_reads_by_bubble_excludes_delete_only_reads`,
    `bubble_site_arm_read_counts_excludes_delete_only_reads`). `cargo clippy
    --all-features --tests` and `cargo fmt --check` clean against this round's
    changes (pre-existing warnings elsewhere in the test suite are unrelated
    and untouched). `PoaConfig`'s defaults are unchanged.

- **`heaviest_path` conflated Match and Delete traversal when scoring the consensus
  spine, letting a node reached mostly by reads *skipping past it* out-compete a
  genuinely Match-confirmed alternative arm at the same fork.** `Node.out_edges`
  stored a single `i32` weight incremented identically by `Match` and `Delete`
  (`increment_or_add_edge` didn't distinguish them), and only node-level
  `coverage`/`delete_count` preserved any Match/Delete split -- as an aggregate
  across *all* of a node's in-edges, losing which specific edge deserved the
  credit. Confirmed on real RFC1 CANVAS data during an architectural audit
  (`design/graph_data_model_rework.md`): nodes reached by one pure-delete in-edge
  and one pure-match in-edge to the same target are not a hypothetical, e.g.
  `(1253, 0 match, 3 delete)` vs. `(2480, 1 match, 0 delete)` both feeding node
  1254. `Node.out_edges` is now `Vec<(usize, EdgeWeight)>` with
  `EdgeWeight { matched: i32, deleted: i32 }` (`Insert`'s founding traversal is
  accounting-identical to a `Match`, so no third bucket is needed); every consumer
  the type change touched got an explicit, documented projection decision rather
  than a blanket behavior change -- `heaviest_path` and the multi-allele bubble
  gates (`find_bubbles`, `find_structural_bubbles`) now score on `matched` only,
  while `GraphStats`, the interior filter (Known Bugs #6-#10), and the public
  diagnostic/visualization APIs (`edge_weights()`, `node_out_edge_info()`,
  `graph_topology()`'s `GraphEdgeInfo.weight`, `bubble_arm_lengths()`) are
  unchanged, still reading total (matched + deleted) weight, deferred to a later
  pass. Matched-only scoring alone was *necessary but not sufficient*: it
  correctly favors the better-evidenced arm at the fork itself, but that local
  advantage can be exactly cancelled at the arms' reconvergence point, because
  reads that Delete through the weak node still genuinely Match whatever follows
  it, and that next edge's own matched weight is real and gets summed in
  regardless -- confirmed by hand on the regression case below, where both arms'
  cumulative scores tied exactly (6 = 6) without an additional fix. Closed by
  also penalizing forward propagation, once per node on the path, by however much
  a node's own `delete_count` exceeds its `coverage` (zero, and so a no-op, for
  the vast majority of ordinary Match-dominant nodes) -- mirroring the interior
  filter's own `coverage > delete_count` philosophy, applied earlier, at
  spine-construction time, rather than only at final-consensus time. Regression:
  `heaviest_path_prefers_matched_over_delete_inflated_arm` (a forced case with
  non-repetitive flanks so semi-global's free boundary gap can't explain the
  length deficit for free: 6 reads with a true deletion allele vs. 4 reads with a
  true SNP allele at the same position). `edge_weight()` (the separate, live
  `align()` M/D-state DP tie-break) was A/B tested both ways against the full
  suite and `bench/validate.py`; neither discriminated, so it stays on total
  weight pending future evidence. Full data-model audit and phased rollout plan
  in `design/graph_data_model_rework.md`; this is Phase 1 of that plan.

- **Multi-allele phasing (`partition_reads_by_bubble`, `phasing_groups`) and
  `BubbleSite.arm_read_counts` inherited the same Match/Delete conflation as
  `heaviest_path` (Phase 1 above), via `PoaGraph.edge_reads`.** `edge_reads`
  recorded read membership for *any* traversal of an edge -- Match, founding
  Insert, or Delete -- with no way to tell them apart, so a read that merely
  deleted through a bubble arm's starting node (skipped it, confirming nothing)
  was indistinguishable from a read that genuinely matched that arm. Both
  `partition_reads_by_bubble` (SNP-bubble allele splitting) and `phasing_groups`
  (structural-bubble cross-compatibility grouping) build each read's per-bubble
  arm signature straight from `edge_reads`, so both inherited the blindness;
  `collect_bubble_sites` reads the same map to populate `BubbleSite.arm_read_counts`.
  `PoaGraph.edge_reads` now records only genuine confirmations (Match, or the
  founding Insert); Delete traversals go into a new, separate
  `edge_delete_reads` map instead -- two parallel maps mirroring the existing
  `coverage`/`delete_count` idiom, rather than one map with a tagged value, so
  every existing call site that reads `edge_reads` for "which reads support this
  arm" gets the corrected meaning with no signature change. Regression:
  `partition_reads_by_bubble_excludes_delete_only_reads` (a forced near-tied
  SNP bubble -- 4 reads confirming each of two arms -- plus one more read
  missing the contested base entirely; before this fix that read was wrongly
  attributed to the arm it deleted through instead of correctly falling
  through as unassigned and folding into the larger group) and
  `bubble_site_arm_read_counts_excludes_delete_only_reads` (same setup,
  confirming `arm_read_counts` no longer counts that read towards either arm).
  Phase 2 of the plan in `design/graph_data_model_rework.md`; Phase 1's
  deferral list (`find_bubbles`/`find_structural_bubbles` weight threshold,
  `GraphStats`, the interior filter, the public diagnostic/visualization APIs)
  is unchanged by this pass.

- **Exact-duplicate divergent arms (a substitution or insertion identical to
  one already created by an earlier read at the same fork) still fragmented
  into separate single-read nodes instead of merging, because nothing indexed
  "the exact sequence of bases this divergence represents" against existing
  arms -- reuse only ever happened by accident, via ordinary DP scoring a
  match against an existing node higher than creating a new one, which does
  not reliably fire in repetitive sequence (see the `LOOKAHEAD_K` finding in
  `design/graph_data_model_rework.md`).** Added a per-fork content-addressable
  arm index, `PoaGraph.fork_arm_index: HashMap<usize, HashMap<Vec<u8>, usize>>`
  (fork node index -> full characterized edit's bases -> existing arm's start
  node index), consulted in `add_to_graph` before creating a new node for a
  `Match` mismatch or an `Insert` run. Keyed by the *complete* edit, not a
  short prefix, specifically to avoid reintroducing the false-merge risk
  `LOOKAHEAD_K`'s own length floor already guards against for the scoring
  path. This is Phase 3 of the plan in `design/graph_data_model_rework.md`;
  it does not replace `LOOKAHEAD_K`/`slide_lock`, and does not attempt
  canonicalization, so fuzzy near-duplicates (the same effective edit
  expressed through incidentally different node chains) are unaffected by
  design, not by oversight.
  - **A naive first version of this idea, tried earlier in the investigation
    that produced the design doc, reused an existing sibling out-edge on a
    matching *first* base alone and caused real hangs (SIGKILL) on real-data
    regression tests** by letting nodes accumulate far more in-degree/fan-in
    than the bounded-depth arm-walking machinery elsewhere in this file
    assumes. This phase's design mitigates that specific failure mode by
    requiring the *entire* edit to match, verified read-only before any
    mutation, before ever reusing a node -- confirmed safe by running the
    same tests that hung before individually, with a bounded timeout, before
    running the full suite.
  - **That mitigation was necessary but not sufficient on its own.** A second,
    different bug surfaced empirically (not assumed fixed by the design):
    reusing an existing node for an `Insert` run redirects the read's
    "current position" to that existing node -- but `align()` computes a
    read's entire traceback in one pass, before `add_to_graph` runs at all,
    and a later `Match`/`Delete` op in that *same* traceback can independently
    name the very node reuse just redirected onto, since the DP had no idea
    reuse would happen. Confirmed concretely on real data (DAB1 SCA37,
    seed=25): node 49 ended up with `in_edges=[17, 49]` and
    `out_edges=[18, 99, 49]`, a literal self-loop, which corrupted
    `topological_order`'s Kahn's-algorithm bookkeeping (in-degree could never
    reach zero) and hung `heaviest_path`'s traceback (an unbounded walk over a
    predecessor-pointer array that develops its own cycle once two different
    nodes end up assigned the same topological rank). Fixed by checking, before
    committing any reuse, whether the candidate chain's nodes are referenced
    again anywhere later in that read's own remaining ops; if so, reuse is
    skipped and a fresh arm is created exactly as before this phase existed.
  - Regression: `period7_content_addressing_reduces_duplicate_forks` (the
    period-7 `GCTAGCT`x10 duplicate-fork audit from the earlier investigation,
    re-measured, not assumed: 166 nodes / 71 singleton(cov=1) nodes / 14
    duplicate-fork positions before this phase, 163 / 67 / 12 after -- a real
    but modest reduction, consistent with this scenario being dominated by
    fuzzy near-duplicates rather than exact ones). All 5 of the bug #6-#10
    `diag_*` regression tests pass individually under a bounded timeout before
    the full suite was run, per the investigation's fail-fast validation order.

- **Diagonal-skip bubble pre-resolve only marked a losing arm's first node as a dead
  end, not the rest of the arm.** When the fast spine-diagonal-skip resolves a fork
  (e.g. a read-supported substitution or short indel creating an alternate node),
  it marks the losing arm's nodes as safe to bypass on future reads. For a 1-node
  arm (a plain substitution) this is the whole arm, so it worked correctly. For a
  2+-node arm (e.g. a 1bp insertion, which creates a two-node detour before
  rejoining the spine), only the arm's first node was ever marked -- every node
  past it was left neither on-spine nor marked, so it silently fell through to a
  real windowed DP computation on every subsequent read, forever, for a branch that
  read will never take. Worse, the arm's reconvergence node kept seeing an
  unresolved incoming edge from the dangling remainder, which failed its own
  fast-path eligibility check and forced every position for the rest of that read
  into full DP too, even when nothing about the rest of the read was ambiguous.
  Fixed by walking the entire losing arm (stopping at the spine or a further
  branch) instead of just its first node, matching what the separate, more
  expensive slide-and-lock bubble resolver already does via `collect_arm_nodes`.
  Regression: `multi_node_minority_arm_fully_marked_as_dead_end`.

- **`PoaConfig::band_width = 0` was not actually unbanded.** Documented as "unbanded
  (full NW over DAG)", but `align()`'s spine-margin computation silently fell back to
  a fixed 50-column window whenever `band_width = 0` and `adaptive_band = false` --
  identical in practice to `band_width = 50`. There was no way to reach the documented
  correctness-baseline behavior through the public API. Found while investigating a
  long-VNTR locus (multi-kb reads, tens-of-bp repeat unit) where the CLI default
  (banded, adaptive on) silently fabricated an extra copy of the repeat unit beyond
  what any input read, abPOA, or SPOA supported on the same reads (a concrete instance
  of the already-documented "silent wrong alignment on narrow band in repetitive
  regions" failure mode: the adaptive formula's band at this read length was just
  wide enough to reach an adjacent, equally-scoring but wrong repeat-phase diagonal,
  while the ~50-wide fallback "unbanded" was actually using was not). `align()`'s
  spine margin now uses the full query length when `band_width = 0` and
  `adaptive_band = false`, matching the documented O(read_len × graph_len) cost and
  recovering the correct repeat count when explicitly requested via
  `--band-width 0 --no-adaptive-band` (the CLI *default* banded behavior on this class
  of locus is unchanged and is tracked separately). `PoaConfig::warn_on_long_unbanded`
  -- documented but never actually checked anywhere -- is now wired up: `align()`
  emits a stderr warning once per read when unbanded DP runs on a read over 1 kb,
  since this is now genuinely expensive rather than silently cheap.
  `consensus_adaptive`'s truncation-retry (which rebuilds with `band_width = 0` when
  the pass-1 consensus looks truncated) gained the same `median_read_len ≤ 5000` cap
  the CLI's own truncation retry already had, so it can no longer trigger an
  unbounded-cost rebuild on long reads now that the retry's fallback is real.
  Regression: `tests/band_width_zero_unbanded.rs` (synthetic fixture; no real patient
  data).

- **The interior filter's fork-context lookup (Known Bugs #6-#10) re-derived
  "is there a fork anywhere back along this node's unbranched run" from scratch
  on every `consensus()` call, via a 64-hop bounded backward walk over the
  heaviest-path spine (`path[search_idx - 1]`, added by bug #8's fix).** This
  worked but was pure re-derivation: the same walk runs again for every node
  in the filtered range, on every call, even though the answer only changes
  when the graph's edge structure changes. Added `Node.nearest_fork:
  Option<(usize, usize)>` (nearest ancestor with 2+ out-edges, and which of
  its out-edges this node descends from), populated incrementally as the
  graph is built rather than walked at consensus time: a brand-new node
  inherits its predecessor's `nearest_fork` unless the predecessor is itself
  a fork, and when a predecessor gains a second out-edge partway through
  the graph's life, `propagate_fork_if_new` walks forward from its existing
  child to backfill `nearest_fork` on every already-existing descendant,
  stopping at a reconvergence point (2+ in-edges), a further fork, or a
  bounded depth (`ARM_MAX_DEPTH`, matching an existing bound already used
  elsewhere in the file, rather than inventing an unproven new one). The
  interior filter's fork lookup itself now reads `Node.nearest_fork` directly;
  its accept/reject judgment (both axes: `coverage > delete_count` when no
  fork is in reach, majority/plurality-of-the-fork's-local-total otherwise)
  is unchanged, confirmed by the full test suite passing byte-for-byte
  identically, not by inspection alone.
  - **The staleness question this design flagged as unproven -- does a
    cache populated incrementally stay correct when a node that was a stable,
    unforked single-successor predecessor for existing descendants later
    becomes a fork -- was real, not hypothetical.** A first version called
    `propagate_fork_if_new` eagerly, inline, at the moment each new edge was
    added. A constructed regression test (several reads establish nodes 5-15
    as a clean single-predecessor chain off node 4, then later reads give
    node 4 a second out-edge) initially failed: propagation correctly walked
    forward and updated node 5, but a downstream node one read's own
    traceback would *later* turn into a genuine reconvergence point was
    updated first, before that same read's remaining ops added the second
    in-edge that should have made it ineligible -- eager propagation had no
    way to see a reconvergence its own read hadn't gotten around to creating
    yet. Fixed by deferring propagation: `add_to_graph` now collects every
    node that gained a new out-edge while processing the *current* read into
    a list, and only calls `propagate_fork_if_new` for each of them once,
    after that read's entire op sequence has been applied -- by which point
    any reconvergence the read itself causes is already reflected in
    `in_edges.len()`. Two regression tests capture this directly:
    `nearest_fork_updates_for_preexisting_descendants_after_late_fork` (the
    non-reconverging case: a permanently divergent one-node arm, confirming
    forward propagation reaches every pre-existing descendant) and
    `nearest_fork_same_read_reconvergence_is_not_stale_updated` (the
    adversarial case: a single read both creates the fork and reconverges
    back onto the original chain a few bases later, confirming the
    reconvergence node and everything past it is correctly left at its old
    value rather than stale-updated) -- both checked against an independent,
    from-scratch backward walk over raw `in_edges`, not against the cache's
    own logic. (White-box tests: `Node.nearest_fork` is a private field, so
    these live in `#[cfg(test)] mod fork_cache_tests` inside `src/graph.rs`
    itself rather than the black-box `src/tests.rs`.)
  - One behavioral-equivalence subtlety needed an explicit guard rather than
    a silent drop: the old walk refused to look further back than the start
    of the current call's boundary-trimmed range (`range_start`), since a
    fork sitting entirely within the trimmed-away leading section shouldn't
    count. `nearest_fork` has no notion of any particular call's boundary
    trim (it's a graph-level property), so the rewired filter reproduces
    the same cutoff explicitly via topological rank comparison
    (`rank_of[pred_idx] >= range_start_rank`) rather than dropping the check
    and letting behavior drift.
  - Phase 4 of the plan in `design/graph_data_model_rework.md`, and the
    highest-risk phase of it by the plan's own assessment -- confirmed by
    this being the first phase where the design's own explicitly-flagged
    "unproven" concern (staleness) turned out to be a real bug, not a false
    alarm.

- **`sv_cag20_out60` (bench validation) failure was almost entirely a test-harness
  artifact, not a `poa-consensus` bug.** `bench/validate.py`'s `build_reference`
  and `simulate_reads` gave every allele in a multi-allele/SV scenario
  byte-identical flanking sequence (both sides sliced from the same offset into
  `_LEFT_FLANK_FULL`/`_RIGHT_FLANK_FULL`), so a read genuinely simulated from one
  allele could cross-map (via `minimap2`) onto a *different* allele's reference
  contig whenever it didn't fully resolve the length-differentiating repeat, and
  `bedpull` would then extract it under the wrong allele label. Confirmed
  concretely: `sv_cag20_out60` (intended: 20 reads of CAG×20 plus a 2-read CAG×60
  outlier) extracted **8** reads tagged `allele_1`, not 2 -- direct inspection
  (read length, and the underlying pbsim read names via `samtools view`) showed
  only 2 were genuinely CAG×60-length; the other 6 were ordinary CAG×20 reads
  that cross-mapped purely because both reference contigs shared identical
  flanks. One of those 6 was short enough that `select_seed(Auto)`'s
  "shortest spanning candidate" heuristic picked it as seed once it entered the
  pool, triggering the same seed-length-sensitivity mechanism described below --
  i.e. the harness bug was what let a spurious short seed into an otherwise
  clean single-allele read set. Fixed by giving each allele index a distinct,
  non-zero rotation offset into the same flank pool (`_rotated_flank`,
  `_FLANK_ROTATE_STEP = 11`, not a multiple of the pool's 32 bp period; verified
  CAG/GAA/CTTTT/AAAAG-free at every offset in range) in both `build_reference`
  and `simulate_reads`'s per-allele pbsim reference. Allele index 0 always gets
  offset 0, so every single-allele scenario (the large majority of the suite)
  and every scenario's first allele are byte-for-byte unaffected. `bench/compare_callers.py`
  imports `build_reference`/`simulate_reads` from `validate.py` directly, so it
  inherits the fix with no separate change needed. Re-running confirmed the
  fix in isolation, with no other change: `sv_cag20_out60` extracts exactly 4
  reads under `allele_1` (matching pbsim's own actual per-allele read count for
  that depth, not the harness's earlier inflated 8), and the CLI's *default*
  auto-seed behavior on the corrected read set alone -- no seed-selection code
  change -- now produces the exact correct CAG×20 consensus (Δ+0, edit=0),
  confirming this specific scenario's failure was **entirely** a harness
  artifact. `bench/validate.py` moved from 17/20 to 18/20 (the remaining 2
  failures, `cag50_d20` and `multi_gaa30_100`, are unrelated -- see the next
  entry for `cag50_d20`'s root cause); `bench/compare_callers.py` stayed 16/16
  with `sv_cag20_out60`/`sv_gaa50_out200` both now showing clean Δ+0 agreement
  across all three callers. Bench-tooling only (`bench/validate.py`); no Rust
  crate changes.

- **Seed-length sensitivity in periodic/homogeneous tandem repeats: an
  auto-selected seed that happens to be atypically short relative to the true
  read population can make the consensus systematically under-call the
  repeat's true length, confirmed on synthetic CAG/GAA data independent of
  band width and independent of the graph-model rework above.** Root-caused via
  ground-truth comparison against pbsim3's own `.maf` simulated-read alignments
  (not fragile substring counting): for a failing `cag50_d20` draw, the actual
  per-read repeat lengths clustered tightly around a true median of 148 bp
  (truth 150 bp), but the consensus recovered only 135 bp (45 of 50 units) --
  13 bp short of *every* read's own true length except the two most
  deletion-heavy outliers. Reproduces identically fully unbanded
  (`band_width = 0, adaptive_band = false`) and at `band_width = 200`, ruling out
  Known Bug #3 (narrow-band diagonal drift) as the mechanism despite the
  superficial resemblance. Mechanism: `select_seed(Auto)` deliberately picks the
  *shortest* spanning candidate (documented, tested behavior -- maximises
  Insert-type bubble visibility); when that read is short relative to the
  population, the majority's extra content has to be inserted somewhere, and in
  a homogeneous repeat any position is an equally valid place to insert it, so
  different reads scatter their inserts across different positions and no
  single insertion position accumulates enough coverage to survive the
  majority/coverage floor on its own. Forcing the same read set to seed on an
  ordinary-length read instead recovered the correct length outright (confirmed
  on multiple independent scenarios and random-seed draws, including through
  the actual CLI binary, not just the library).
  - **Investigated as a candidate absolute "is this consensus suspicious"
    signal, and rejected, with the negative results being the useful part of
    this entry:** `GraphStats::single_support_fraction`, `bubble_count`,
    `edge_weight_gini`, `mean_column_entropy`, chosen-seed length relative to
    the read population's median/longest read, and a new aggregate "off-spine
    fork mass" measure (built for this investigation: for every spine fork,
    what fraction of its total traffic goes to a losing arm, aggregated across
    the graph) were all tested across every confirmed-failing case and 21
    confirmed-passing scenarios spanning CAG/GAA at multiple lengths and
    depths, ONT R9/R10, and HiFi. None separated cleanly: e.g.
    `single_support_fraction` sits in the same 0.29-0.46 range for essentially
    every periodic-repeat scenario tested, whether anything is wrong or not,
    and a passing scenario (`cag50_d30_s42`) had a *higher*
    `single_support_fraction` (0.458) than any of the three confirmed-failing
    cases. This confirms the *existing* `consensus_adaptive` trigger on this
    same statistic (`single_support_fraction > 0.3`) was never actually
    detecting "this consensus is wrong" -- only "this locus looks
    noisy/repetitive," which is true of nearly all of them.
  - **What did discriminate, empirically: comparing several candidate
    consensuses built from the *same* reads against each other**, scored by a
    new function, `analysis::consensus_fit(consensus_seq, reads, config)` --
    builds a throwaway graph on `consensus_seq` alone, aligns every read
    against it, and reports the mean per-read total of Insert ops (content the
    candidate doesn't explain) plus Delete ops (candidate content the read
    doesn't confirm), normalised by length. Lower is a better fit. This is
    explicitly a *relative* scorer across candidates from one read population,
    not an absolute per-scenario threshold -- validated as such, not as a
    standalone "is this wrong" check. Two variants were tried: `mean_insert_frac`
    alone (more aggressive: fixed all 3 known failures in testing, but
    regressed 2 of 21 passing scenarios, one catastrophically -- a
    longest-length seed candidate scored deceptively well because a long,
    error-inflated seed reduces how much *other* reads need to insert, without
    reducing how much they need to *delete*, which this variant doesn't
    penalise); and the symmetric Insert+Delete version actually shipped, which
    never chose a worse candidate than pass-1 for any of 21 tested passing
    scenarios in this investigation.
  - **Wired into `consensus_adaptive`'s existing `single_support_fraction > 0.3`
    branch** (`src/lib.rs`), replacing its previous unconditional remedy
    (tighten `min_coverage_fraction` to >= 0.6 and rebuild on the same seed --
    kept as one candidate, since it does help in some cases, e.g.
    `cag60_d20_s42` in this investigation) with a scored comparison across four
    candidates: pass-1 unchanged, the pre-existing tightened-coverage rebuild,
    a rebuild re-seeded on a read near the population's median length
    (`median_length_read_index`, new), and a `ConsensusMode::MajorityFrequency`
    rebuild (same seed) -- keeping whichever scores best by `consensus_fit`.
    Extracted the scoring logic into a private `seed_sensitivity_retry` helper
    shared with a new standalone public entry point,
    `consensus_fit_scored(reads, seed_idx, config)`, for callers that want just
    this retry without `consensus_adaptive`'s other passes (multi-allele bubble
    detection, truncation retry, semi-global fallback) -- this crate's own CLI
    (`src/main.rs`) now calls it for its single-allele path instead of plain
    `consensus()`, since the CLI already had a documented reason to avoid
    `consensus_adaptive` itself (multi-allele bubbles firing unexpectedly on a
    het SNP in a caller that expects single-allele output).
  - **Validated end to end, including through the compiled CLI binary, not just
    the library:** full `cargo test --release --features cli` (one existing
    test, `adaptive_noisy_tightens_coverage`, needed its `action` assertion
    relaxed from asserting `NoisyTighten` unconditionally to accepting either
    `NoisyTighten` or `PassThrough` -- confirmed this is the new design working
    as intended, not a regression: the scattered single-base-error scenario it
    constructs was already correctly resolved by pass-1's ordinary majority
    vote with no tightening needed, and `consensus_fit` now correctly detects
    that and returns pass-1 unchanged instead of needlessly rebuilding; the
    *sequence* assertion, which is what the test actually cares about, is
    unchanged and still passes). `cargo clippy --all-features --tests` and
    `cargo fmt --check` clean, matching the pre-change warning baseline
    exactly. `bench/validate.py`: 18/20, unchanged from the harness-fix-only
    baseline -- `cag60_d20_s42` (an ad hoc seed-sweep case from this
    investigation, not a tracked scenario) went from Δ-8 to Δ+0 through the
    real CLI; `cag50_d20` (tracked, named scenario) is unaffected either way,
    since its own best-scoring candidate ties with pass-1; one ad hoc sweep
    case (`cag50_d20_s2`) that was *already* failing before this change went
    from Δ-7 to Δ-10 -- confirmed not a regression of a previously-passing
    scenario, but reported plainly since it is a real, acknowledged limit of
    this fix, not swept under the rug. `bench/compare_callers.py`: 16/16,
    unchanged. Not a complete fix for the underlying seed-length sensitivity --
    see `consensus_adaptive`'s doc comment for the full, honest account of what
    this does and does not solve.

- **`bench/validate.py`'s `extract_allele_by_unit` was itself fabricating most
  of the apparent severity of `cag50_d20` and (per the entry above)
  `cag50_d20_s2` -- both look far closer to correct once measured properly.**
  Digging into why `consensus_fit` scored the (supposedly under-called) pass-1
  consensus as tied with the true answer for `cag50_d20`, and why the retry's
  chosen candidate for `cag50_d20_s2` scored *better* than pass-1 despite the
  benchmark reporting it as *more* wrong (Δ-7 -> Δ-10), turned up the same root
  cause behind both: `extract_allele_by_unit`'s cluster-gap tolerance
  (`2*len(unit)`, 6 bp for `CAG`/`GAA`) is too tight. A single interior
  substitution error (e.g. one `CAG` -> `CAA`) breaks exact-match scanning for
  a stretch wide enough to push the gap between surviving hits past 6 bp with
  no actual indel present at all, splitting one genuine repeat run into two
  clusters; the function then keeps only the larger cluster and silently
  discards the rest as "not really part of the repeat."
  - **Confirmed independently of this function, not just by adjusting its own
    tolerance and re-checking itself:** built the ground-truth window directly
    from the reference-construction constants (100 bp of real flank on each
    side + the true repeat) and compared it to the actual consensus output via
    plain Levenshtein distance (no clustering, no unit-counting at all).
    `cag50_d20`'s pass-1 consensus (reported as Δ-5 units/edit=15) is only
    **4 edits** from the 350 bp ground-truth window -- essentially exact.
    `cag50_d20_s2`'s retry-chosen candidate (reported as Δ-10/edit=30, i.e.
    "worse than pass-1's Δ-7") is only **8 edits** from ground truth, against
    pass-1's own 28 -- the retry was a real, substantial improvement, not a
    regression; the previous entry's "Δ-7 -> Δ-10" framing is superseded by
    this direct measurement. Also re-ran `consensus_fit` on these consensuses
    via the actual crate function (`poa_consensus::analysis::consensus_fit`),
    not a re-implementation: its scores for `cag50_d20` (pass-1 0.050883 vs.
    ground truth 0.050795) and its ranking of all four `cag50_d20_s2`
    candidates matched the Levenshtein-based ordering exactly in both cases --
    `consensus_fit` was never the problem; it was scoring correctly against a
    ground truth that the benchmark's own reporting function was
    misrepresenting.
  - **Fixed by widening the tolerance from `2*len(unit)` to `4*len(unit)`.**
    Verified safe before applying to the real file: patched a copy of the
    function in memory and re-evaluated all 16 tracked single-allele scenarios
    plus 9 ad hoc seed-sweep cases from this investigation (25 total) --
    zero pass/fail status changes anywhere. Applied to the real file and
    re-confirmed the same 25/25 zero-flip result directly against it (not the
    copy). `bench/validate.py` stays 18/20 and `bench/compare_callers.py` stays
    16/16 -- **this changes reported severity, not pass/fail status**:
    `cag50_d20` now reports Δ-2/edit=2 (was Δ-5/edit=15) and still fails the
    strict ±1-unit tolerance, `cag50_d20_s2` now reports Δ-5/edit=7 (was
    Δ-10/edit=30) and still fails, and `cag20_d20_r9` (already excused via the
    read-limit adequacy check) improves from Δ-7/edit=21 to Δ-2/edit=2 with no
    status change. Both real residuals (Δ-2 and Δ-5) are now small enough that
    whether they reflect a genuine remaining algorithmic gap or further
    benchmark-measurement noise is an open question for future investigation,
    not something this fix resolves -- it only removes a confound that was
    making both failures look roughly 5-10x more severe than the consensus
    output actually was. Bench-tooling only (`bench/validate.py`); no Rust
    crate changes; `cargo test`/`clippy`/`fmt` are unaffected (confirmed, as
    expected for a Python-only change).

### Breaking

- **`BubbleSite.arm_read_counts` now counts only reads that genuinely confirmed
  an arm (Match, or the founding Insert), not reads that merely deleted through
  its entry node.** The field's type is unchanged (`Vec<u32>`), but its value
  can differ from earlier versions for any bubble where a read skipped past an
  arm's starting node without confirming it -- such a read no longer inflates
  that arm's count. Counting only genuine confirmations is the corrected
  behavior (see the Match/Delete conflation entry under Fixed above); callers
  relying on the old "any traversal" semantics, e.g. to sum `arm_read_counts`
  as a proxy for total depth at a bubble, should account for the reduced
  counts.

- **`AdaptiveAction` gained two new variants, `AlternateSeedRetry` and
  `MajorityFrequencyRetry`, and `single_support_fraction > 0.3` no longer
  guarantees `NoisyTighten`.** The enum is not `#[non_exhaustive]`, so an
  exhaustive `match` on `AdaptiveAction` in downstream code will fail to
  compile until the two new arms are added. Semantically: before this change,
  `consensus_adaptive` returned `NoisyTighten` unconditionally whenever
  `single_support_fraction > 0.3` fired; it now builds several candidate
  remedies (see the seed-length-sensitivity entry under Fixed above), scores
  them against the actual reads, and returns whichever action produced the
  best-scoring candidate -- which may be `PassThrough` (pass-1 was already
  best), the pre-existing `NoisyTighten`, or one of the two new variants.
  Callers that specifically branched on `action == NoisyTighten` to mean
  "the trigger fired" should check for all four outcomes (or none, and just
  use `consensuses[0]`) instead.

## [0.2.1] - 2026-07-08

### Fixed

- **Lookahead bubble-arm scoring biased toward the longer arm** -- the lookahead resolver
  (used at bubble entries with 2+ out-edges to pick a winning arm before windowed DP) gated
  on whether the *longest* arm could provide `LOOKAHEAD_K` bases, then scored each arm only
  over its own available length. A short arm (e.g. a 1-node "no insertion" shortcut) could
  therefore only ever score up to `1 x match_score`, while a long arm (e.g. a 19-node
  insertion detour created by one divergent read) scored over the full `LOOKAHEAD_K` bases --
  guaranteeing the long arm won the `LOOKAHEAD_MARGIN` comparison regardless of which arm
  actually matched the read. Once one outlier read (with a genuinely different repeat-unit
  count) created such a bubble in a periodic homopolymer/STR run, every subsequent read
  hitting that bubble was routed onto the same over-scored long arm, silently fragmenting
  read support across the two arms and producing a consensus with several small deletions
  scattered through the repeat -- with no error and only an easy-to-miss low-coverage
  warning. The result was highly seed-index-dependent: some seeds avoided the outlier
  entirely (correct output), others produced mild (1 bp) to severe (10+ bp) base loss.
  `slide_lock` already guarded against this (`min_arm_len < LOOKAHEAD_K -> bail`); the
  lookahead gate now uses the same *minimum*-arm-length check, so it only fires when every
  arm can be scored over the same number of bases. Confirmed on DAB1 SCA37 (ATTTT/ATTTC)
  hap2 HiFi reads.
- **DP tie-break now prefers the heavier-weighted edge** -- in the windowed-DP fallback (used
  when the lookahead/slide-lock arm resolvers decline), ties between candidate predecessors
  with identical scores were broken by `in_edges` iteration order (edge creation order), not
  by any property of the read data. Ties now prefer the predecessor edge with higher existing
  weight, so genuinely ambiguous alignments reinforce whichever path already has more read
  support instead of fragmenting further. This is a secondary hardening alongside the
  lookahead fix above; it did not by itself resolve the DAB1 case but closes the same class of
  order-dependent tie-breaking.
- **CLI default seed selection** -- `poa-consensus` (no `--seed`) now uses
  `select_seed(&reads, &SeedSelection::Auto)` (terminal k-mer anchor heuristic, already used
  internally by the library) instead of a local median-read-length heuristic. The median
  heuristic had no particular basis for correctness and, before the lookahead fix above, made
  the CLI's default behaviour on the DAB1 case worse than several other seed choices.
- **Interior filter could fabricate sequence not present in any read** -- `consensus()`'s
  interior filter dropped any node with `coverage < min_cov`, where `min_cov` is derived from
  the *global* read count. A node that `heaviest_path` correctly picked as the winning arm of
  a local bubble can still have Match coverage below that global threshold, purely because
  other reads diverged at an *earlier* bubble (so fewer than the full read count even reach
  this point) -- the node is a clear local majority but a global minority. Dropping it spliced
  its flanking nodes directly together, producing a junction no read ever made. This is acute
  in periodic/tandem-repeat regions: removing one interior base of a `CTT` repeat shifts the
  phase (`"CTTCTTC"` with the middle `C` removed reads as `"CTTTTC"`, a run length the repeat
  unit can never actually produce). The interior filter now also keeps a node when it is the
  locally dominant arm at its own bubble -- its Match coverage clears a majority of its
  predecessor's total out-edge weight -- gated on the predecessor actually having 2+ out-edges
  (a real fork), so a plain low-coverage unbranched tail is still governed by the global
  threshold alone. Confirmed on DMD `CTT`-repeat HiFi reads (HG02968, chrX:31,284,557-613);
  regression test: `diag_dmd_ctt_repeat_interior_filter_global_threshold_bias`.
  **Known residual** (not fixed by the above): when the *majority* of reads reach a node via
  Delete rather than Match -- i.e. they genuinely skip the base, not just diverge onto another
  arm -- dropping that node can still occasionally fabricate a junction, because each deleting
  read's own alignment still passes through the node in sequence and no single edge directly
  connects its predecessor to its successor bypassing it. An initial attempt to also require a
  same-length bypass edge before dropping any node fixed this case but broke the (correct)
  removal of genuine multi-node minority insertions elsewhere (regressed
  `diag_dab1_sca37_attttc_lookahead_arm_length_bias`), so it was not kept. Tracked by the
  `#[ignore]`d test `diag_dmd_ctt_majority_delete_residual`.
- **Local-dominance rescue (above) could itself fabricate sequence from a fragmented bubble's
  noise** -- the rescue gates only on *relative* local majority (a node's Match coverage
  clearing a majority of its predecessor's local out-edge weight total), not on the *absolute*
  size of that local population. In a region with many small, closely-spaced phase-registration
  bubbles (confirmed on RFC1 `AAAAG` pentanucleotide repeat, HG002 real data), the read pool
  fragments rapidly bubble-by-bubble; by the time one specific fork is reached, only a handful
  of reads may still be on that exact sub-path. A 2-of-3 split there clears the local majority
  bar trivially even though it's statistical noise from attrition, not a real minority allele
  the global threshold unfairly suppresses -- reported as a spurious single-base `G` interrupt
  inserted after unit 30 of a 116-unit `AAAAG` run, with no read support at that position (2 of
  10 reads had an `AAAAGG` anywhere in their whole sequence, both far from unit 30). Fixed by
  additionally requiring the bubble's local population (`pred`'s total out-edge weight) to clear
  the *global* `min_cov` before trusting a local majority within it. Regression test:
  `diag_rfc1_aaaag_spurious_g_interrupt_local_rescue_noise` (sweeps all 10 possible seeds).
- **Local-dominance rescue only checked the *immediate* predecessor for a fork, missing the
  case one hop downstream of an already-rescued node** -- confirmed on RFC1 `AAAAG`, HG002
  real data, hap2: a read set with 4 full-length spanning reads plus 2 short partial reads that
  only cover the first ~55-60bp. Once the 2 partial reads' alignments ended (a legitimate
  termination under semi-global mode, not a bubble), `n_reads` stayed pinned at 6 for the
  *global* `min_cov`, even though only 4 reads structurally reach deeper positions. A node
  immediately downstream of an already-rescued fork -- whose own predecessor has only one
  out-edge, i.e. isn't itself a fork -- always fell through to the global check and got dropped
  purely because some reads elsewhere in the graph ended earlier, even when it was the clear
  majority (3 of the 4 still-active reads). Fixed by walking backward (bounded to 64 hops) to
  find the nearest real fork rather than requiring the immediate predecessor to be one, so a
  rescued node's own single-successor chain inherits its fork's population instead of being
  judged against the global read count. Regression test:
  `diag_rfc1_leading_interrupt_partial_read_population_accounting` (checks the region confirmed
  fully explained by population accounting via a 4-full-reads-only control).
  Also surfaced a pre-existing, unrelated issue while testing this data: using one of the two
  short partial reads as the seed produces `TruncationRetry` with `recovered: false` and an
  empty consensus (confirmed present before this session's changes too); not investigated
  further here.
- **Backward-fork-search rescue (above) still missed the case where no fork exists anywhere
  nearby** -- confirmed on the same RFC1 `AAAAG` hap2 read set: one read individually `Delete`d
  a single `G` out of an otherwise fully unbranched, unanimous run of matches. No fork is ever
  created for this, because Match and Delete share the same edge -- only the node's own
  `delete_count` records that one read skipped it. The node's own Match coverage (3 of the 4
  reads that structurally reach this deep) fell below the *global* `min_cov` (6, sized for all
  10 reads including the 2 short ones that never reach this point) purely because no fork
  existed to rescue it via the previous fix. Dropping that `G` merged the two flanking 4-`A`
  groups either side of it into one 8-`A` run, reported as an `(A)4(AAAAG)1(AAAG)1` interrupt
  with zero read support at that position -- this was the residual reported as "the same class
  of fabrication further into the repeat" in the previous entry; it turned out to be a third
  instance of the same population-accounting bug, not the periodic-alignment-ambiguity issue it
  was first mistaken for. Fixed by falling back, when no fork is found in the backward search,
  to the node's own `coverage + delete_count` as the local population -- every read reaching an
  unforked node either matches or deletes it, so that sum is exact -- gated on the sum itself
  clearing the *global* `min_cov`, mirroring the fork branch's own gate.
  **This gate is load-bearing, not optional**: an earlier version of this fix used a
  forward-recovery scan instead (does coverage climb back above `min_cov` further along?) and
  regressed two existing tests -- `long_repeat_length_majority_wins` and
  `edge_extreme_length_variation_majority_wins` in `tests/sv_analysis.rs` -- both genuine
  trailing extensions (a minority of reads simply longer than the rest, with nothing downstream
  to reconverge with, so nothing should be rescued). The forward scan could be fooled by a
  single stray noisy read overlapping the tail; comparing the population itself against the
  global threshold cannot, because a genuine trailing minority's own population is, by
  definition, still a minority. Regression test:
  `diag_rfc1_leading_interrupt_delete_driven_forkless_gap`.
  With this fix, the RFC1 hap2 read set is now clean end-to-end except for a single, much
  smaller residual (`gap=2` at one position) that *is* the genuine periodic-alignment-ambiguity
  class (Known Bug #3/#4): a minority read's private insertion absorbing Match support from
  other reads' independent re-alignment, since any `A` matches any `A` node identically. Not
  addressed by this fix; flanking-anchor pre-processing remains the intended long-term fix.
- **`bench/validate.py` used removed CLI flags** -- `run_consensus()` passed `--adaptive-band`
  and `--semi-global`, both removed when those became the CLI defaults (see "CLI changes"
  below). Every `validate.py` invocation was silently producing `FAILED (no consensus)` for
  every scenario. Updated to pass `--band-width 50` (single-allele) / `--band-width 0`
  (multi-allele, pure adaptive) with no other band/mode flags.

### Changed

- **Interior filter simplified: two evidence axes instead of one threshold ladder.** Prompted
  by `bench/compare_callers.py` surfacing a low-depth gap (`cag20_d05_r10`, depth 5) where
  abPOA and SPOA both matched truth exactly and poa-consensus did not. Reading abPOA's
  (`abpoa_heaviest_bundling`/`abpoa_most_frequent` in `abpoa_output.c`) and SPOA's
  (`Graph::TraverseHeaviestBundle`/`GenerateConsensus` in `graph.cpp`) actual source showed
  neither applies any absolute population floor to a DP-selected node at all -- our own
  `increment_or_add_edge` call site confirmed Match and Delete increment the *same* edge
  weight identically, so `heaviest_path`'s own DP score already conflates "many reads reached
  this point" with "many reads support this exact base," and the interior filter had grown
  into an increasingly complex ladder (global floor, then backward fork search, then a
  self-total fallback with its own floor) trying to patch that conflation back apart after
  the fact. Replaced with two independent, uniform tests, applied per spine node after a
  boundary trim unchanged from before (an absolute `coverage >= min_cov` floor is still the
  only way to distinguish a genuine interior dip from a trailing/leading minority extension,
  since the latter has `delete_count == 0` too -- the rest of the population simply never
  reaches that far under semi-global alignment, so a Match-vs-Delete comparison alone can't
  tell the two apart):
  1. **Match vs Delete** (does this exact base have more support than "skip it"): applies
     when no fork (2+ out-edges) exists anywhere back along the node's unbranched run, all the
     way to the trim boundary. In that case coverage and delete_count together are the exact,
     complete count of reads that reached this point, so plain `coverage > delete_count` is
     sufficient -- no population floor needed. This is the axis majority-delete fabrication
     (Known Bugs #6, #9) and a false rescue of one (traced to DAB1 SCA37 seed=0 mid-implementation:
     a node with cov=3, del=23 was wrongly kept by an edge-weight-based plurality check, because
     edge weight counts Match and Delete traversals together) both come down to.
  2. **This base vs a different base** (a real fork somewhere back along the run): needs the
     fork's own coverage to clear a *majority* of the fork's local total, gated on that local
     total itself clearing the global floor (Known Bug #7's fragmented-fork guard, unchanged) --
     *unless* both the fork node itself and the candidate arm have `delete_count == 0`, in which
     case a plurality (the arm whose entry weight is at least as large as every sibling's) is
     trusted instead. A bare plurality relaxation (no delete_count check at all) was tried and
     reverted twice first: it wrongly rescued a genuine 5-vs-5 near-tie in real RFC1 AAAAG data
     (fabricating an unsupported extra base between two branches of the same fork -- but that
     fork's *own* node had `cov=4, del=6`, i.e. the arrival at the fork was itself an unresolved
     majority-delete decision one level up), and it wrongly rescued 10-of-14 reads' insertion arm
     against a clean 4-unit-majority CAT tandem-duplication test (the rejection at the fork's
     first node wasn't propagating to the rest of that arm's unbranched chain, fixed separately
     by having the backward search below track which of the fork's direct children the whole arm
     descends from, so one verdict covers the entire arm, not just its first node). Gating the
     plurality relaxation on `delete_count == 0` on *both* ends -- confirmed by checking the
     RFC1 case's actual node data -- distinguishes "a fork with a genuinely undisputed arrival,
     splitting cleanly into several different-but-real bases" (`cag20_d05_r10`: fork cov=7,
     del=0; candidate cov=3, del=0 -- safe for plurality) from "a fork whose own arrival is still
     contested" (RFC1: fork cov=4, del=6 -- not safe, needs strict majority regardless of the
     split below it).
  Net effect on `bench/compare_callers.py`'s 16-scenario catalogue: 16/16 all three callers now
  agree (up from 15/16) -- `cag20_d05_r10` genuinely fixed, `sv_cag20_out60` unaffected. All 184
  lib tests plus the full `tests/` suite still pass with zero failures, despite fixing DAB1 by a
  cleaner path, finding + fixing the tandem-duplication regression, and fixing the
  originally-targeted low-depth gap for real.

### Added

- **`bench/compare_callers.py`** -- runs the same pbsim3 -> minimap2 -> bedpull pipeline as
  `validate.py`, but feeds the extracted reads to `poa-consensus`, abPOA (`pyabpoa`), and SPOA
  (`pyspoa`) side by side, all in semi-global mode, and scores each against the known truth
  allele. Restricted to single-allele scenarios: neither external tool has an equivalent to
  this crate's bubble-based multi-allele splitting, so a fair three-way comparison needs all
  three solving the same one-read-set-in-one-consensus-out problem. Setup: `pip install
  pyabpoa pyspoa`. Scoring went through two broken extraction-based attempts before landing on
  `fitting_edit_distance()` (semi-global alignment of the truth allele against the full
  consensus, no extraction step at all) -- see `TODO.md` Deferred for why both extraction
  approaches failed on this benchmark's synthetic data. The comparison surfaced the low-depth
  gap resolved by the interior filter changes above; full run against the existing scenario
  catalogue now agrees 16/16 across all three callers.

---

## [0.2.0] - 2026-06-22

### Fixed

- **Network diagram arm layout** (`plot` feature) -- multi-node insertion arms (e.g. a 2-base
  insertion) are now grouped into chains before y-slot assignment, so all nodes in the same
  arm land on the same side of the spine. Nodes within a chain are spaced one spine-step apart
  horizontally so they do not overlap. Previously the slot counter incremented per node,
  causing the second node of a 2-node arm to flip to the opposite side.
- **Shifted minimizer anchor bug in `anchor_refine_spine`** -- k-mers that span a
  repeat/flanking boundary can appear at a different read position in shorter reads (one fewer
  repeat unit before the anchor sequence). This produced anchors whose upper bound `ahi` fell
  below the expected alignment diagonal `j_center`, constraining `j_hi` below the reachable
  range and causing catastrophic alignment failure: the affected read produced M=1, I=99
  (roughly 100 spurious insert nodes) instead of the correct M=100, D=3. `anchor_refine_spine`
  now discards any anchor where `ahi < j_center` rather than applying it, leaving the band
  at its unrefined default. Confirmed on SCA4 ZFHX3 hap2 HiFi reads; regression test:
  `diag_sca4_zfhx3_hap2_alignment_ops`.

### Added

- **`graph_network_svg`** (`plot` feature) -- new visualisation of POA graph topology as a
  directed network. Spine nodes are laid out horizontally; arm (insert) nodes sit above or
  below their attachment points; edges between spine and arm nodes are drawn as curves. Labels
  are rendered inside each node. An optional read slice overlays that read's alignment path in
  a distinct colour, making it easy to inspect which nodes a specific read traverses. Run the
  included example with `cargo run --example network_plot --features plot`.
- **`graph_network_svg_labeled`** (`plot` feature) -- same as `graph_network_svg` but labels
  each edge with its read count (edge weight). Spine edges use kuva's `with_edge_label`;
  arm edge weights are injected directly into the SVG after rendering, above or below each
  arm node circle. Useful for illustrating the heaviest-path DP, boundary trim thresholds,
  and alignment debugging.
- **`Consensus::read_indices: Vec<usize>`** -- indices into the original `reads` slice that
  contributed to this consensus. Populated by `consensus_multi` and `PoaGraph::consensus_multi`;
  empty for single-allele outputs (empty means "all reads contributed"). Callers can use these
  to assign per-read rows to the correct haplotype, pull reads for visualisation, or compute
  per-allele statistics without re-running alignment.
- **`AdaptiveAction` enum and `AdaptiveResult` struct** -- `consensus_adaptive` now returns
  `Result<AdaptiveResult, PoaError>` instead of `Result<Vec<Consensus>, PoaError>`.
  `AdaptiveResult` carries `consensuses: Vec<Consensus>` and `action: AdaptiveAction`.
  `AdaptiveAction` variants: `PassThrough`, `MultiAllele`, `NoisyTighten`,
  `TruncationRetry { recovered: bool }`, `SemiGlobalFallback`. Callers that previously
  used the return value as a `Vec` directly must now access `.consensuses`.

### Breaking

- **`consensus_adaptive` return type** changed from `Result<Vec<Consensus>, PoaError>` to
  `Result<AdaptiveResult, PoaError>`. Replace `consensus_adaptive(...)?` with
  `consensus_adaptive(...)?.consensuses` at call sites that only need the consensus vec.

---

## [0.1.2] - 2026-06-01

### Changed

- **`PoaConfig::default()` now uses `AlignmentMode::SemiGlobal`** (was `Global`). Semi-global
  alignment (free terminal gaps on the graph) is correct for extracted STR reads, which start
  and end at slightly different positions depending on where each read clips into the flanking
  sequence. Global alignment silently produces wrong boundary behaviour when the seed has more
  leading identical bases than most reads: shorter reads place their deletions inside the
  homopolymer run (not at the start), leaving the extra node on the only spine path where
  boundary trim cannot reach it. Confirmed on FRDA FXN HiFi data.
- **`PoaConfig::default()` now uses `adaptive_band: true, band_width: 50`** (was `false, 0`).
  The previous unbanded default was impractical for reads above ~1 kb and silently produced
  wrong unit counts in long repetitive alleles when called with the recommended adaptive + 50
  floor settings. The new default matches the recommended single-allele STR configuration.
  The 50 bp floor prevents unit loss that the raw adaptive formula (w ≈ 15 for 500 bp reads)
  would allow. Callers that need unbanded DP (e.g. correctness tests, short reads where memory
  is not a concern) can set `band_width: 0, adaptive_band: false` explicitly.

### CLI changes

- **`--semi-global` removed**; semi-global is now the default. Use `--global` to opt into
  global alignment (useful when reads are guaranteed to span the full locus from identical
  start/end positions).
- **`--adaptive-band` removed**; adaptive band is now the default. Use `--no-adaptive-band`
  to disable.
- **`--band-width` default changed from `0` to `50`**. Together with adaptive band, the
  effective band is `max(50, 10 + 0.01 × read_len)`. Pass `--no-adaptive-band --band-width 0`
  for fully unbanded alignment.

---

## [0.1.1] - 2026-06-01

### Added

- **Truncation detection** -- `diagnose()` now sets `ConsensusWarnings::truncation_suspected`
  (`Option<TruncationWarning>`) when the consensus length is less than
  `DiagnoseConfig::truncation_ratio_threshold` (default `0.60`) times the median input read
  length. `TruncationWarning` carries `consensus_len`, `median_read_len`, and `ratio` for
  downstream inspection.
- **`GraphStats::median_input_read_len`** -- new field; set automatically by `consensus()` and
  `stats()` from the reads stored in the graph. Required by the truncation detector.
- **`DiagnoseConfig::truncation_ratio_threshold`** -- new field (default `0.60`). Set to `0.0`
  to disable truncation detection entirely.
- **CLI unbanded retry** -- the single-allele path in `poa-consensus` now checks `diagnose()`
  after the initial banded run. When truncation is suspected and the median read length is
  ≤ 5,000 bp, it automatically retries with `band_width = 0` (unbanded). For longer reads,
  where an unbounded DP matrix is impractical, it emits a warning suggesting `--band-width 0`
  unless `--quiet` is set.
- **`bench/hg002_real_data.py`** -- HG002 real-data validation script covering eight disease
  STR loci (HTT, FMR1, DMPK, RFC1, ATXN1, ATXN2, ATXN3, and ATXN3 with corrected
  coordinates) plus three non-repetitive controls (GAPDH, TP53 exon 7, ACTB). Features:
  CIGAR-based read clipping, HP-tag splitting for haplotagged BAMs, per-locus
  `poa_extra_flags` override with `merge_poa_flags()` deduplication, bimodal read-length
  estimator, and a truth table with tolerance-aware evaluation.

### Fixed

- **RFC1 HiFi consensus** (bug #4 workaround) -- the AAAAG 5-mer alignment degeneracy that
  caused silent graph fragmentation on HiFi reads is now handled by the truncation-detection
  unbanded retry. Using `pad=20` (20 bp unique flanking context) makes the fragmented
  consensus fall below the 0.60 ratio threshold; the unbounded retry recovers ~106× AAAAG
  (HG002 truth ~115×, within ±10 tolerance). This is a workaround; flanking-anchor
  pre-processing remains the intended long-term fix.
- **GAPDH single-allele regression** -- the CLI single-allele path no longer routes through
  `consensus_adaptive`, which fired the multi-allele early-return on heterozygous SNPs in
  non-repetitive loci (producing a 22 bp truncated output). The path now calls `consensus()`
  directly followed by an independent truncation check.

---

## [0.1.0] - 2026-06-01

Initial release.

### Added

- Banded affine-gap POA: tracking band, diagonal skip, bubble-range j-window, stale-spine
  adaptive recompute, approaching-edge detector, smart 3-pass retry, lookahead resolve,
  exact-size slide-lock DP for committed bubble arms.
- Single-allele (`consensus`), multi-allele (`consensus_multi`), and adaptive two-pass
  (`consensus_adaptive`) functional wrappers over the stateful `PoaGraph` API.
- `SeedSelection::Auto` -- terminal k-mer anchor heuristic; selects the shortest spanning
  read, falling back gracefully when no spanning read exists.
- `AlignmentMode::Global` and `AlignmentMode::SemiGlobal` (free terminal gaps for partial
  reads).
- `ConsensusMode::HeaviestPath` and `ConsensusMode::MajorityFrequency`.
- Structural bubble phasing via union-find over pairwise arm-traversal agreement; falls back
  to single-best SNP bubble when no structural bubbles are found.
- `diagnose()` / `ConsensusWarnings` / `DiagnoseConfig` -- actionable diagnostics: low depth,
  coverage gaps, interior low support, competing structural alleles.
- Analysis helpers: `consensus_confidence`, `has_competing_allele`, `should_call_multiallele`,
  `allele_fractions`, `count_credible_interval`, `max_achievable_accuracy`,
  `low_coverage_regions`, `min_coverage`.
- `extract_flanked_region` -- flanking-anchor pre-processing for rotation-anchored repeat
  extraction.
- `bridged_consensus` -- joins two partial-read consensuses with a `GapKind::Unknown` gap.
- Orientation utilities: `reverse_complement`, `orient_to_seed`, `auto_orient`.
- `GraphStats` -- O(V+E) summary: bubble count, node coverage, edge weight Gini coefficient,
  single-support fraction, longest bubble span, and more.
- CLI (`--features cli`): FASTA/FASTQ input via noodles, clap argument parsing, `--multi`,
  `--semi-global`, `--adaptive-band`, `--band-width`, `--quiet`, seed selection.
- Plot helpers (`--features plot`): coverage SVG, graph-stats bar chart, edge-weight
  histogram, node-coverage histogram, alignment-density heatmap, band-with-reads overlay.
- 193 tests; 20/20 synthetic validation scenarios passing (two via adequacy signals).

[Unreleased]: https://github.com/Psy-Fer/poa-consensus/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/Psy-Fer/poa-consensus/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/Psy-Fer/poa-consensus/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/Psy-Fer/poa-consensus/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/Psy-Fer/poa-consensus/releases/tag/v0.1.0
