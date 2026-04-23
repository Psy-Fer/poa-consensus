use std::fmt;

#[derive(Debug)]
pub enum PoaError {
    EmptyInput,
    InsufficientDepth { got: usize, min: usize },
    SeedOutOfBounds { index: usize, len: usize },
    BandTooNarrow { band_width: usize, read_len: usize },
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
                band_width,
                read_len,
            } => {
                write!(
                    f,
                    "band width {band_width} is too narrow for read length {read_len}; \
                     alignment reached the band edge and may be incorrect"
                )
            }
        }
    }
}

impl std::error::Error for PoaError {}
