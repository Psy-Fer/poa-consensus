# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [0.1.3] - unreleased

### Fixed

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
  arm edges get a `TextAnnotation` positioned near the arm node. Useful for illustrating the
  heaviest-path DP, boundary trim thresholds, and alignment debugging.
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

[Unreleased]: https://github.com/Psy-Fer/poa-consensus/compare/v0.1.3...HEAD
[0.1.3]: https://github.com/Psy-Fer/poa-consensus/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/Psy-Fer/poa-consensus/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/Psy-Fer/poa-consensus/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/Psy-Fer/poa-consensus/releases/tag/v0.1.0
