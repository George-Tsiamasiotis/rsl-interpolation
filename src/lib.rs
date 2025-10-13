//! A re-write of [`GSL's Interpolation`] in Rust.
//!
//! [`GSL's Interpolation`]: https://www.gnu.org/software/gsl/doc/html/interp.html
//!
//! # 1D Interpolation Types
//!
//! + [Linear]
//! + [Cubic]
//! + [CubicPeriodic]
//! + [Akima]
//! + [AkimaPeriodic]
//! + [Steffen]
//! + [DynInterpType] - Dynamic Interpolation Type
//!
//! # 2D Interpolation Types
//!
//! + [Bilinear]
//! + [Bicubic]
//!
//! # Higher level Interface
//!
//! + [Spline]
//! + [Spline2d]
//! + [DynSpline]
//! + [DynSpline2d]
//!
#![allow(rustdoc::broken_intra_doc_links)]
#![doc = include_str!("../TODO.md")]

mod accel;
mod error;
mod interp;
mod interp2d;
mod types;

mod spline;
mod spline2d;

pub use accel::Accelerator;

pub use error::*;
pub use interp::{InterpType, Interpolation};
pub use interp2d::{Interp2dType, Interpolation2d, z_get, z_idx, z_set};

pub use spline::{DynSpline, Spline, make_spline};
pub use spline2d::{DynSpline2d, Spline2d, make_spline2d};

pub use types::*;

#[cfg(test)]
mod tests;

/// Trait for supported data types.
pub trait Num: num::Float + num_traits::NumAssignOps + Default + Sync + Send {}

impl Num for f64 {}
impl Num for f32 {}
