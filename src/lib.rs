//! A re-write of [`GSL's Interpolation`] in Rust.
//!
//! [`GSL's Interpolation`]: https://www.gnu.org/software/gsl/doc/html/interp.html

mod accel;
mod error;
mod interp_types;

#[cfg(test)]
mod gsl_tests;

pub use accel::Accelerator;
pub use error::InterpError;
pub use interp_types::InterpolationType;
