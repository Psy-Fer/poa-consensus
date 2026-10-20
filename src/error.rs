use std::fmt;

/// Errors returned by the consensus entry points and [`select_seed`](crate::select_seed).
#[derive(Debug)]
pub enum PoaError {
    /// No reads were provided.
    EmptyInput,
    /// Fewer reads than [`PoaConfig::min_reads`](crate::PoaConfig::min_reads)
    /// (in multi-allele mode: too few non-empty reads in total, or in one
    /// allele group).
    InsufficientDepth {
        /// Reads available.
        got: usize,
        /// The `min_reads` floor that was not met.
        min: usize,
    },
    /// The seed index is not a valid index into the reads.
    SeedOutOfBounds {
        /// The index passed.
        index: usize,
        /// Number of reads.
        len: usize,
    },
    /// `SeedSelection::Auto` found non-overlapping left-only and right-only
    /// read groups with no spanning read.  Use `bridged_consensus` instead:
    /// the left group seeds the left half, the right group seeds the right half,
    /// and `bridged_consensus` concatenates them with a `GapKind::Unknown` gap.
    NoSpanningReads {
        /// Reads that cover only the left side.
        left_depth: usize,
        /// Reads that cover only the right side.
        right_depth: usize,
    },
}

impl fmt::Display for PoaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PoaError::EmptyInput => {
                write!(f, "no reads provided")
            }
            PoaError::InsufficientDepth { got, min } => {
                write!(
                    f,
                    "insufficient depth: got {got} reads, need at least {min}"
                )
            }
            PoaError::SeedOutOfBounds { index, len } => {
                write!(f, "seed index {index} is out of bounds for {len} reads")
            }
            PoaError::NoSpanningReads {
                left_depth,
                right_depth,
            } => {
                write!(
                    f,
                    "no spanning reads found: {left_depth} left-only and {right_depth} \
                     right-only reads detected — use bridged_consensus to assemble each \
                     side separately"
                )
            }
        }
    }
}

impl std::error::Error for PoaError {}
