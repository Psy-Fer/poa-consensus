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
    /// builds with it `true`; single-allele paths leave it `false`). Stateful
    /// callers ([`PoaGraph::new`]) who intend to call
    /// [`PoaGraph::consensus_multi`] should set it `true` so alignment is built
    /// in multi-allele mode. Default: `false` (single-allele).
    ///
    /// [`consensus_multi`]: crate::consensus_multi
    /// [`PoaGraph::new`]: crate::PoaGraph::new
    /// [`PoaGraph::consensus_multi`]: crate::PoaGraph::consensus_multi
    pub multi_allele: bool,
    /// Add an abPOA-style log-odds **read-support term** to the alignment DP.
    ///
    /// In a homogeneous tandem repeat, matching a base to a phase-shifted repeat
    /// node scores identically (`+match_score`) to matching the correct node, so
    /// the aligner can silently fold a read onto the wrong repeat unit and
    /// fabricate a phantom base at the flank/repeat boundary (Known Bug #3; see
    /// `design/flank_fabrication_bug3.md`). This is a *scoring degeneracy*, and
    /// no band width or centering fixes it. When enabled, every DP transition
    /// across an edge `p -> n` is biased by
    /// `round(log(edge_matched_weight / total_out_matched_weight(p)))` -- the log
    /// of the fraction of reads that took that edge, `<= 0`, clamped at `-20` --
    /// so a read prefers the well-supported (modal) diagonal on a tie. Mirrors
    /// abPOA's `--inc-path-score` / `-G` option.
    ///
    /// Only takes effect on the **single-allele** path (`multi_allele == false`):
    /// biasing alignment toward the heaviest existing arm would corrupt
    /// multi-allele bubble phasing. Off by default (like abPOA's `-G`); the CLI
    /// enables it for its single-allele consensus path. Costs a `ln()` per DP
    /// transition, so it is opt-in.
    pub path_score_bias: bool,
}

impl Default for PoaConfig {
    fn default() -> Self {
        PoaConfig {
            band_width: 50,
            adaptive_band: true,
            adaptive_band_b: 10,
            adaptive_band_f: 0.01,
            match_score: 1,
            mismatch_score: -1,
            gap_open: -2,
            gap_extend: -1,
            min_coverage_fraction: 0.0,
            min_allele_freq: 0.2,
            min_reads: 1,
            alignment_mode: AlignmentMode::SemiGlobal,
            consensus_mode: ConsensusMode::HeaviestPath,
            warn_on_long_unbanded: true,
            phasing_bubble_min_span: 10,
            multi_allele: false,
            path_score_bias: false,
        }
    }
}
