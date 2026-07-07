# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

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
