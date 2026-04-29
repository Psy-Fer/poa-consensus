use std::fmt;

#[derive(Debug)]
pub enum PoaError {
    EmptyInput,
    InsufficientDepth { got: usize, min: usize },
    SeedOutOfBounds { index: usize, len: usize },
    BandTooNarrow { configured: usize, required: usize },
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
            PoaError::BandTooNarrow { configured, required } => {
                write!(
                    f,
                    "band width {configured} too narrow; estimated {required} required — \
                     retry with a wider band or enable adaptive_band"
                )
            }
        }
    }
}

impl std::error::Error for PoaError {}
