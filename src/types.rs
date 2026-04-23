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
    /// Per-base coverage along the consensus path.
    pub coverage: Vec<u32>,
    pub graph_stats: GraphStats,
}
