/// Orientation of a read relative to a reference read, as returned by
/// [`orient_to_seed`](crate::orient_to_seed).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    /// Same strand as the reference (also returned when there is no k-mer evidence).
    Forward,
    /// Opposite strand: reverse-complement the read before POA.
    Reverse,
}

/// Whether the size of a [`CoverageGap`] can be estimated from the data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapKind {
    /// A spanning read crossed this region.  `sequence[start..end]` are
    /// seed-only bases; `end - start` is a minimum size estimate.
    Spanning,
    /// No read crossed this region.  `start == end` (no bases in the output
    /// represent the gap); true size is completely unknown.
    Unknown,
}

/// A region in the consensus where coverage dropped to seed-only level, or a
/// structural break between two non-overlapping read groups.
///
/// ```text
/// kind = Spanning  →  |===[ seed-only bases ]===|   (≥N bp, N = size())
/// kind = Unknown   →  |===|?????|===|              (at least left + right bp)
/// ```
///
/// Typical rendering:
/// ```text
/// format!("{}(gap:{}{}bp){}", &seq[..gap.start],
///     if gap.kind == GapKind::Unknown { "unknown, ≥" } else { "≥" },
///     gap.size(), &seq[gap.end..])
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoverageGap {
    /// First position in the consensus sequence (0-indexed, inclusive) that
    /// belongs to the gap.  For `Unknown` gaps `start == end`.
    pub start: usize,
    /// First position after the gap (0-indexed, exclusive).
    /// For `Unknown` gaps `start == end`.
    pub end: usize,
    /// Whether a spanning read provided a minimum size estimate.
    pub kind: GapKind,
}

impl CoverageGap {
    /// Number of seed-only bases in the gap (minimum size estimate).
    /// Returns 0 for `Unknown` gaps where `start == end`.
    pub fn size(&self) -> usize {
        self.end - self.start
    }

    /// Minimum number of bases in this gap, if known.
    /// Returns `None` for `GapKind::Unknown`.
    pub fn min_size(&self) -> Option<usize> {
        match self.kind {
            GapKind::Spanning => Some(self.end - self.start),
            GapKind::Unknown => None,
        }
    }
}

/// One bubble site in the POA graph: a node on the heaviest path that has two
/// or more outgoing arms each confirmed by at least
/// `max(1, ceil(n_reads * min_allele_freq))` reads.
///
/// Arms with fewer reads (but above the frequency threshold) are competing
/// alleles or sequencing artefacts; callers use `arm_read_counts` and
/// `arm_sequences` to decide which.  Sites are found on the heaviest path even
/// when [`ConsensusMode`](crate::ConsensusMode) emits the majority-frequency
/// sequence instead.
///
/// `arm_sequences` is empty for arms longer than the internal cap (`ARM_SEQUENCE_CAP`
/// = 256 bp); use `arm_read_counts` to assess support for large structural variants.
///
/// # Example
/// ```text
/// // 10 reads: 7 took the CAG×3 arm, 3 took the CAG×4 arm.
/// site.arm_read_counts  == [7, 3]
/// site.arm_sequences    == [b"CAGCAGCAG", b"CAGCAGCAGCAG"]
/// site.is_structural    == true   // length-changing
///
/// // The 3-read arm exceeds min_allele_freq=0.25 → consider re-running
/// // with consensus_multi to separate the two alleles.
/// ```
#[derive(Debug, Clone)]
pub struct BubbleSite {
    /// Position in the consensus sequence (0-indexed) of the entry node just
    /// before the bubble branches.  Use this to locate the site in the output.
    pub consensus_pos: usize,
    /// Number of reads that genuinely confirmed each arm (reached the arm's
    /// first node directly from the branching node with no intervening
    /// deletion), in arm order.
    ///
    /// A read that merely *deleted through* the branch point (skipped nodes
    /// without confirming any base there) is **not** counted here, so a
    /// length variant's deletion bypass does not show up as a spurious arm.
    /// Because the consensus path follows total edge weight (which does
    /// include such reads), the highest-count arm is usually, but not always,
    /// the one on the consensus.
    pub arm_read_counts: Vec<u32>,
    /// Sequence of each arm (from the first arm node up to but not including the
    /// exit/reconvergence node).  Empty `Vec` when the arm exceeds 256 bp or is
    /// a direct edge to the exit (0-length deletion arm).
    pub arm_sequences: Vec<Vec<u8>>,
    /// `true` when at least one arm spans `≥ phasing_bubble_min_span` bases.
    /// Structural bubbles are length-changing variants (indels, repeat-count
    /// differences); non-structural bubbles are SNPs or short MNPs.
    pub is_structural: bool,
}

/// Per-graph statistics computed in a single O(V+E) pass.
#[derive(Debug, Clone, Default)]
pub struct GraphStats {
    /// Number of nodes in the graph.
    pub node_count: usize,
    /// Number of edges in the graph.
    pub edge_count: usize,
    /// Number of fork nodes anywhere in the graph with 2+ outgoing arms each
    /// confirmed by at least `max(1, ceil(n_reads * min_allele_freq))` reads.
    pub bubble_count: usize,
    /// Largest minority-arm read count over those forks (the read count of the
    /// second-best arm at the most balanced bubble), not a nesting depth.
    /// 0 when there is no qualifying bubble.
    pub max_bubble_depth: usize,
    /// Mean coverage across all nodes.
    pub coverage_mean: f64,
    /// Coverage variance across all nodes.
    pub coverage_variance: f64,
    /// Gini coefficient of edge weights (0 = uniform, 1 = all weight on one edge).
    pub edge_weight_gini: f64,
    /// Fraction of nodes supported by exactly one read.
    pub single_support_fraction: f64,
    /// Mean per-column Shannon entropy (bits).
    pub mean_column_entropy: f64,
    /// Span in bases of the longest arm across all bubbles that meet the
    /// `min_allele_freq` threshold.  0 when no qualifying bubble exists.
    /// Arms longer than 4096 nodes are capped at 4096.
    pub longest_bubble_span: usize,
    /// Median length (in bases) of all reads that built this graph, including
    /// the seed.  0 when no reads have been added.  Used by [`diagnose`](crate::diagnose) to
    /// detect consensus truncation: a consensus much shorter than the median
    /// input read is a signal that banded DP converged to the wrong diagonal.
    pub median_input_read_len: usize,
}

/// Output of a consensus extraction pass.
#[derive(Debug, Clone)]
pub struct Consensus {
    /// The consensus sequence.
    pub sequence: Vec<u8>,
    /// Per-base count of reads occupying that position (node coverage on the
    /// heaviest path, or column occupancy for a majority-frequency sequence).
    pub coverage: Vec<u32>,
    /// Per-base support for the step into each base.  On the heaviest path:
    /// the incoming edge weight (the first base uses its node coverage).  When
    /// the majority-frequency sequence is emitted (forced, or chosen by
    /// best-fit): `min(coverage[i-1], coverage[i])`, with `coverage[0]` for
    /// the first base.
    pub path_weights: Vec<i32>,
    /// Total reads used to build this consensus (including the seed).
    pub n_reads: usize,
    /// Statistics of the graph this consensus was extracted from.
    pub graph_stats: GraphStats,
    /// Coverage gaps detected in this consensus.  Empty when reads overlap
    /// throughout.  Each gap is either `Spanning` (seed-only bases, minimum
    /// size known) or `Unknown` (no spanning read; size completely unknown,
    /// as produced by [`bridged_consensus`](crate::bridged_consensus)).
    pub gaps: Vec<CoverageGap>,
    /// Bubble sites on the consensus path where two or more arms each have
    /// read support above `min_allele_freq`.  Empty when all forks are
    /// single-arm (no competing allele above threshold).
    ///
    /// Each site corresponds to one branching node on the heaviest path.  Inspect `arm_read_counts` and `arm_sequences` to decide
    /// whether to re-run with [`consensus_multi`](crate::consensus_multi).
    pub bubble_sites: Vec<BubbleSite>,
    /// Indices of the reads that contributed to this consensus, as indices into
    /// the input `reads` slice you passed (regardless of `seed_idx`).  Populated
    /// by the [`consensus_multi`](crate::consensus_multi) free function; **empty
    /// for single-allele outputs** ([`consensus`](crate::consensus)).
    ///
    /// An empty `Vec` means "all reads contributed" — no phasing was performed.
    /// A non-empty `Vec` identifies exactly which reads belong to this allele,
    /// enabling the caller to assign per-read rows, pull reads for visualisation,
    /// or compute per-allele statistics without re-running alignment.
    pub read_indices: Vec<usize>,
    /// Per-base support quality, Phred-scaled and capped at 60: how decisively
    /// the reads back each consensus base over its strongest alternative.
    ///
    /// Each read votes from its route through the graph: following the
    /// consensus into a base agrees with it; any other route (a substitution,
    /// a skip, an extra insertion) counts against the bases it bypasses. Reads
    /// taking the same alternative route are pooled, and each base is scored
    /// only against its strongest alternative, so scattered sequencing errors
    /// do not add up but a competing call does. With `a` reads for the base and
    /// `m` for its strongest alternative, the value is
    /// `-10 log10 P(Binomial(a + m + 1, 1/2) > a)`: about Q33 for 10 agreeing
    /// reads with no alternative, Q60 (the cap) from about 20, Q6 for a 6-4
    /// split and Q3 for a tie or no spanning reads.
    ///
    /// It measures agreement, not truth: an error every read shares (such as
    /// a systematic homopolymer undercall) still scores high. Use it to mark
    /// solid and weak stretches of a call, e.g. with
    /// [`low_support_regions`](crate::analysis::low_support_regions) (Q10 is
    /// a good cut).
    pub support_quality: Vec<u8>,
}

impl Consensus {
    /// Per-base fraction of reads supporting each consensus base, in [0.0, 1.0].
    ///
    /// The numerator is `path_weights[i]` (the number of reads that traversed
    /// the incoming edge at this position).  The denominator is `n_reads`.
    /// For partial-read inputs, positions covered only by the spanning seed
    /// will have weight 1 and a low fraction.
    pub fn weight_fraction(&self) -> Vec<f32> {
        if self.n_reads == 0 {
            return vec![0.0; self.path_weights.len()];
        }
        self.path_weights
            .iter()
            .map(|&w| (w as f32 / self.n_reads as f32).clamp(0.0, 1.0))
            .collect()
    }
}
