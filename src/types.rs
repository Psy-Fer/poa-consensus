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
