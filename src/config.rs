#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentMode {
    Global,
    SemiGlobal,
}

/// How the consensus sequence is extracted from the built graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConsensusMode {
    /// Default. "Best-fit": build both a heaviest-path and a majority-frequency
    /// consensus and keep whichever the reads better support (lower mean
    /// per-read insert+delete on realignment). The two win on different inputs
    /// — heaviest on clean/short, majority on high-error length-variable repeats
    /// — so picking per call gets the better of both. Never worse than plain
    /// heaviest path.
    HeaviestPath,
    /// Force the majority-frequency (MSA-column) consensus regardless of fit:
    /// each column emits its plurality base, counting read deletions explicitly.
    /// Best when column majority is trusted outright, e.g. high-depth amplicons.
    MajorityFrequency,
}

#[derive(Debug, Clone)]
pub struct PoaConfig {
    /// 0 = unbanded (full NW over DAG). Warn if reads exceed ~1 kb and this is 0.
    pub band_width: usize,
    /// Enable adaptive band width using the abPOA formula: w = b + f*L.
    pub adaptive_band: bool,
    /// Base component of the adaptive band formula (abPOA default: 10).
    pub adaptive_band_b: usize,
    /// Length-proportional component of the adaptive band formula (abPOA default: 0.01).
    pub adaptive_band_f: f32,
    pub match_score: i32,
    pub mismatch_score: i32,
    /// One-time penalty when a gap opens. Negative.
    pub gap_open: i32,
    /// Per-base penalty inside a gap. Negative.
    pub gap_extend: i32,
    /// Fraction of reads that must cover a node for it to appear in consensus.
    pub min_coverage_fraction: f64,
    /// Minimum fraction of reads supporting an allele for bubble detection.
    pub min_allele_freq: f64,
    /// Minimum number of reads required to build a consensus.
    pub min_reads: usize,
    pub alignment_mode: AlignmentMode,
    pub consensus_mode: ConsensusMode,
    /// Emit a warning to stderr when reads exceed ~1 kb with band_width = 0.
    pub warn_on_long_unbanded: bool,
    /// Minimum arm span (nodes) for a bubble to be considered a structural variant
    /// for cross-bubble phasing. Bubbles below this threshold (SNPs, short indels)
    /// use the existing single-bubble partitioning. Default 10.
    pub phasing_bubble_min_span: usize,
    /// Whether this graph is being built for **multi-allele** consensus.
    ///
    /// Controls two mode-dependent behaviours whose correct setting differs
    /// between single- and multi-allele consensus:
    /// - the O(1) diagonal-skip fast path in alignment (kept in multi-allele
    ///   mode to lock reads onto their own length-allele track; disabled in
    ///   single-allele mode where its greedy forward-match over-calls periodic
    ///   repeats by matching reads through phantom units), and
    /// - the whole-graph unbanded rebuild on band-retry (needed only for
    ///   multi-allele bubble-structure consistency).
    ///
    /// The functional wrappers set this automatically ([`consensus_multi`]
    /// builds with it `true`; single-allele paths leave it `false`), so most
    /// callers never touch it. Default: `false` (single-allele).
    ///
    /// [`consensus_multi`]: crate::consensus_multi
    pub multi_allele: bool,
}

impl Default for PoaConfig {
    fn default() -> Self {
        PoaConfig {
            band_width: 50,
            adaptive_band: true,
            adaptive_band_b: 10,
            adaptive_band_f: 0.01,
            // abPOA-style scoring: gaps and mismatches are harsh relative to
            // match (+2), so the aligner does not open cheap spurious gaps in
            // homopolymer/periodic runs. The old +1/-1/-2/-1 made gaps too cheap
            // (a gap-open cost only 2 matches), scattering homopolymer alignments
            // at high error and over-calling repeats; abPOA-like scoring roughly
            // halves that error with no regression on clean inputs (validated on
            // the robustness matrix + 3-way comparison, 2026-07-28).
            match_score: 2,
            mismatch_score: -4,
            gap_open: -4,
            gap_extend: -3,
            min_coverage_fraction: 0.0,
            min_allele_freq: 0.2,
            min_reads: 1,
            alignment_mode: AlignmentMode::SemiGlobal,
            consensus_mode: ConsensusMode::HeaviestPath,
            warn_on_long_unbanded: true,
            phasing_bubble_min_span: 10,
            multi_allele: false,
        }
    }
}
