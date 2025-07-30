//! A re-write of [`GSL's Interpolation`] in Rust.
//!
//! [`GSL's Interpolation`]: https://www.gnu.org/software/gsl/doc/html/interp.html

mod accel;
mod error;
mod interp;

mod types;

pub use error::*;
pub use interp::Interpolation;

pub use accel::Accelerator;
pub use types::*;

#[cfg(test)]
mod gsl_tests;
