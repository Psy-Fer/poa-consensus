use std::fmt;

#[derive(Debug)]
pub enum PoaError {
    EmptyInput,
    InsufficientDepth {
        got: usize,
        min: usize,
    },
    SeedOutOfBounds {
        index: usize,
        len: usize,
    },
    BandTooNarrow {
        configured: usize,
        required: usize,
    },
    /// `SeedSelection::Auto` found non-overlapping left-only and right-only
    /// read groups with no spanning read.  Use `bridged_consensus` instead:
    /// the left group seeds the left half, the right group seeds the right half,
    /// and `bridged_consensus` concatenates them with a `GapKind::Unknown` gap.
    NoSpanningReads {
        left_depth: usize,
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
            PoaError::BandTooNarrow {
                configured,
                required,
            } => {
                write!(
                    f,
                    "band width {configured} too narrow; estimated {required} required — \
                     retry with a wider band or enable adaptive_band"
                )
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
