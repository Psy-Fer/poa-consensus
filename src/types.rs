#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
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
