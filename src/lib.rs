pub mod alexander_table;
pub mod error;
pub mod geometry;
pub mod hull;
pub mod io;
pub mod kmt;
pub mod knotsize;
pub mod knottype;
pub mod point;
pub mod polynomial;

// Re-export primary public API
pub use alexander_table::AlexanderTable;
pub use error::{KnotError, Result};
pub use knotsize::{find_knot_core, KnotCoreResult, KnotSizeOptions};
pub use knottype::{get_knottype, get_knottype_ring, KnotOptions};
pub use point::Point3;
