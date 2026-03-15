use thiserror::Error;

#[derive(Debug, Error)]
pub enum KnotError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("failed to parse polynomial: {0}")]
    PolynomialParse(String),

    #[error("failed to parse data: {0}")]
    DataParse(String),

    #[error("convex hull computation failed")]
    HullFailed,

    #[error("knot type not found for computed polynomial: {0}")]
    NotFound(String),

    #[error("empty chain provided")]
    EmptyChain,
}

pub type Result<T> = std::result::Result<T, KnotError>;
