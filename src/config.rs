#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentMode {
    Global,
    SemiGlobal,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConsensusMode {
    HeaviestPath,
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
}

impl Default for PoaConfig {
    fn default() -> Self {
        PoaConfig {
            band_width: 0,
            adaptive_band: false,
            adaptive_band_b: 10,
            adaptive_band_f: 0.01,
            match_score: 1,
            mismatch_score: -1,
            gap_open: -2,
            gap_extend: -1,
            min_coverage_fraction: 0.0,
            min_allele_freq: 0.2,
            min_reads: 1,
            alignment_mode: AlignmentMode::Global,
            consensus_mode: ConsensusMode::HeaviestPath,
            warn_on_long_unbanded: true,
        }
    }
}
