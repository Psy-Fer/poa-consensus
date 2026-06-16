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
    /// Number of reads that traversed each arm, in arm order.  The arm with the
    /// highest count is the one chosen for the consensus sequence.
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
