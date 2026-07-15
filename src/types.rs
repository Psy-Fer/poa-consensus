#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
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

/// One bubble site in the POA graph: a node on the consensus path that has two
/// or more outgoing arms each supported by at least `min_allele_freq` of reads.
///
/// The heaviest arm determines the consensus sequence.  Arms with fewer reads
/// (but above the frequency threshold) are competing alleles or sequencing
/// artefacts; callers use `arm_read_counts` and `arm_sequences` to decide which.
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
    /// Number of reads that genuinely confirmed each arm (Match, or the
    /// founding Insert that created it), in arm order.  The arm with the
    /// highest count is the one chosen for the consensus sequence.
    ///
    /// A read that merely *deleted through* an arm's starting node -- skipped
    /// it without confirming any base there -- is **not** counted here, even
    /// though it structurally traversed the same edge. Before this crate's
    /// `EdgeWeight`/`edge_reads` Match/Delete split (see
    /// `design/graph_data_model_rework.md`), this field conflated the two, so
    /// a delete-heavy skip-through could inflate an arm's apparent support.
    /// Counting only genuine confirmations is a behavior change from earlier
    /// versions, even though the field's type is unchanged (`Vec<u32>`).
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

/// A single node in the POA graph, as returned by [`PoaGraph::graph_topology`].
#[derive(Debug, Clone)]
pub struct GraphNodeInfo {
    /// Internal node index (not stable across calls to `add_read`).
    pub node_idx: usize,
    /// DNA base at this node (`b'A'`, `b'C'`, `b'G'`, or `b'T'`).
    pub base: u8,
    /// Number of reads with a Match op at this node.
    pub coverage: u32,
    /// Number of reads that passed through via Delete (traversed without consuming a base).
    pub delete_count: u32,
    /// Topological position of this node (0 = first node).
    pub topo_rank: usize,
}

/// A directed edge in the POA graph, as returned by [`PoaGraph::graph_topology`].
#[derive(Debug, Clone)]
pub struct GraphEdgeInfo {
    /// Topological rank of the source node.
    pub from_rank: usize,
    /// Topological rank of the target node.
    pub to_rank: usize,
    /// Number of reads that traversed this edge.
    pub weight: i32,
}

/// A snapshot of the POA graph topology for visualization and inspection.
///
/// Returned by [`PoaGraph::graph_topology`].  Indices in `edges` refer to
/// entries in `nodes` by `topo_rank`.
#[derive(Debug, Clone)]
pub struct GraphTopology {
    /// Nodes in topological order.
    pub nodes: Vec<GraphNodeInfo>,
    /// All directed edges.
    pub edges: Vec<GraphEdgeInfo>,
    /// Topological ranks of nodes on the heaviest-path (consensus) spine.
    pub spine_ranks: Vec<usize>,
}

/// Per-graph statistics computed in a single O(V+E) pass.
#[derive(Debug, Clone, Default)]
pub struct GraphStats {
    pub node_count: usize,
    pub edge_count: usize,
    /// Number of heterozygous nodes (nodes with 2+ successors each above min_allele_freq).
    pub bubble_count: usize,
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
    /// the seed.  0 when no reads have been added.  Used by [`diagnose`] to
    /// detect consensus truncation: a consensus much shorter than the median
    /// input read is a signal that banded DP converged to the wrong diagonal.
    pub median_input_read_len: usize,
}

/// Output of a consensus extraction pass.
#[derive(Debug, Clone)]
pub struct Consensus {
    pub sequence: Vec<u8>,
    /// Per-base match-read count along the consensus path.
    pub coverage: Vec<u32>,
    /// Incoming edge weight on the heaviest path for each base.
    /// For the first base (no incoming path edge) the node's coverage is used.
    /// For `ConsensusMode::MajorityFrequency` every position uses node coverage.
    pub path_weights: Vec<i32>,
    /// Total reads used to build this consensus (seed + all `add_read` calls).
    pub n_reads: usize,
    pub graph_stats: GraphStats,
    /// Coverage gaps detected in this consensus.  Empty when reads overlap
    /// throughout.  Each gap is either `Spanning` (seed-only bases, minimum
    /// size known) or `Unknown` (no spanning read; size completely unknown,
    /// as produced by [`bridged_consensus`]).
    pub gaps: Vec<CoverageGap>,
    /// Bubble sites on the consensus path where two or more arms each have
    /// read support above `min_allele_freq`.  Empty when all forks are
    /// single-arm (no competing allele above threshold).
    ///
    /// Each site corresponds to one branching node on the heaviest-path
    /// consensus.  Inspect `arm_read_counts` and `arm_sequences` to decide
    /// whether to re-run with [`PoaGraph::consensus_multi`].
    pub bubble_sites: Vec<BubbleSite>,
    /// Indices (into the original `reads` slice) of the reads that contributed
    /// to this consensus.  Populated by [`PoaGraph::consensus_multi`] and the
    /// [`consensus_multi`] free function; **empty for single-allele outputs**
    /// ([`PoaGraph::consensus`], [`consensus`]).
    ///
    /// An empty `Vec` means "all reads contributed" — no phasing was performed.
    /// A non-empty `Vec` identifies exactly which reads belong to this allele,
    /// enabling the caller to assign per-read rows, pull reads for visualisation,
    /// or compute per-allele statistics without re-running alignment.
    pub read_indices: Vec<usize>,
}

/// The action taken by [`consensus_adaptive`] on its second pass.
///
/// Returned inside [`AdaptiveResult`] so callers can distinguish a clean
/// pass-through from a corrected result without re-running [`diagnose`].
///
/// [`consensus_adaptive`]: crate::consensus_adaptive
/// [`diagnose`]: crate::diagnose
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AdaptiveAction {
    /// Pass-1 result was returned unchanged. Either no second-pass trigger
    /// fired at all, or the singleton-support trigger fired but the pass-1
    /// consensus itself scored best among the candidate remedies compared
    /// (see [`NoisyTighten`](AdaptiveAction::NoisyTighten)) -- i.e. the
    /// pass-1 result was already trustworthy despite the trigger.
    PassThrough,
    /// Competing bubble(s) above `min_allele_freq` triggered multi-allele
    /// partitioning.  The `consensuses` vec has one element per detected allele.
    MultiAllele,
    /// Pass-1 consensus was below the truncation-ratio threshold; alignment was
    /// retried with unbanded DP (`band_width = 0`).
    TruncationRetry {
        /// `true` if the unbanded retry produced a consensus at or above the
        /// truncation threshold.  `false` means the output is still short:
        /// the sequence may be genuinely short or the locus may need manual
        /// review.
        recovered: bool,
    },
    /// High singleton-support fraction triggered a comparison of several
    /// candidate remedies, scored empirically against the actual read
    /// population (see [`consensus_fit`](crate::analysis::consensus_fit)),
    /// and the pass-1 consensus rebuilt with `min_coverage_fraction` raised
    /// to ≥ 0.6 scored best.
    ///
    /// Before the seed-sensitivity retry (below) was added, this variant
    /// was returned unconditionally whenever the trigger fired; it is now
    /// one of several scored candidates, so the trigger firing no longer
    /// guarantees this specific action -- see
    /// [`AlternateSeedRetry`](AdaptiveAction::AlternateSeedRetry) and
    /// [`MajorityFrequencyRetry`](AdaptiveAction::MajorityFrequencyRetry)
    /// for the other possible outcomes, and
    /// [`PassThrough`](AdaptiveAction::PassThrough) for the case where the
    /// pass-1 consensus itself was already the best-scoring candidate.
    NoisyTighten,
    /// High singleton-support fraction triggered the same empirical
    /// candidate comparison as [`NoisyTighten`](AdaptiveAction::NoisyTighten),
    /// and re-seeding on a different read (one near the read population's
    /// median length, rather than the auto-selected seed) scored best.
    ///
    /// Confirmed root cause this targets: an auto-selected seed that is
    /// atypically short relative to the true read population can cause
    /// `heaviest_path`/the interior filter to systematically under-call a
    /// periodic/homogeneous repeat's length, because the extra content the
    /// majority of reads carry (relative to the short seed) gets inserted at
    /// ambiguous, scattered positions and no single insertion accumulates
    /// enough coverage to survive on its own. Re-seeding on a more
    /// representative-length read avoids the problem outright rather than
    /// trying to recover from it after the fact.
    AlternateSeedRetry,
    /// High singleton-support fraction triggered the same empirical
    /// candidate comparison as [`NoisyTighten`](AdaptiveAction::NoisyTighten),
    /// and rebuilding with [`ConsensusMode::MajorityFrequency`] (same seed)
    /// scored best.
    ///
    /// [`ConsensusMode::MajorityFrequency`]: crate::ConsensusMode::MajorityFrequency
    MajorityFrequencyRetry,
    /// High coverage coefficient of variation in `Global` mode triggered a
    /// rebuild with [`AlignmentMode::SemiGlobal`].
    ///
    /// [`AlignmentMode::SemiGlobal`]: crate::AlignmentMode
    SemiGlobalFallback,
}

/// Return value of [`consensus_adaptive`].
///
/// [`consensus_adaptive`]: crate::consensus_adaptive
#[derive(Debug, Clone)]
pub struct AdaptiveResult {
    /// Assembled consensus sequences.  One element for single-allele outcomes;
    /// two or more for multi-allele.
    pub consensuses: Vec<Consensus>,
    /// Which second-pass action (if any) was taken.
    pub action: AdaptiveAction,
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
