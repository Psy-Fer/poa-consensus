# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Fixed

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
