//! A re-write of [`GSL's Interpolation`] in Rust.
//!
//! ## Notes
//!
//! 1. `rsl-interpolation` requires LAPACK FFI, so you must use **just one** of the corresponding
//! [`ndarray_linalg features`]:
//!
//! ```text
//! [dependencies]
//! rsl-interpolation = { version = "0.1.19", features = ["ndarray-linalg/openblas-system"]}
//! ```
//!
//! 2. In 2d Interpolation, the `za` array must be defined in **column-major (Fortran)** style.
//! This is done to comply with GSL's interface.
//!
//!
//! # 1D Interpolation Types
//!
//! + [`LinearInterpolator`]
//! + [`CubicInterpolator`]
//! + [`CubicPeriodicInterpolator`]
//! + [`AkimaInterpolator`]
//! + [`AkimaPeriodicInterpolator`]
//! + [`SteffenInterpolator`]
//!
//! # 2D Interpolation Types
//!
//! + [`BilinearInterpolator`]
//! + [`BicubicInterpolator`]
//!
//! # Lower level interface
//!
//! + Interpolation type must be known at compilation time.
//! + The [`Interpolation`]/[`Interpolation2d`] objects do not store the arrays they were
//! constructed with; they must be passed as arguments at each evaluation.
//!
//! ## 1D Interpolation
//!
//! ```
//! # use rsl_interpolation::*;
//! # fn main() -> Result<(), InterpolationError>{
//! let xa = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
//! let ya = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
//! let interp = CubicInterpolator::build(&xa, &ya)?;
//!
//! let x = 3.9;
//! let acc = &mut Accelerator::new();
//! let y = interp.eval(&xa, &ya, x, acc)?;
//! let dy = interp.eval_deriv(&xa, &ya, x, acc)?;
//! let dy2 = interp.eval_deriv2(&xa, &ya, x, acc)?;
//! let int = interp.eval_integ(&xa, &ya, 0.4, 3.1, acc)?;
//! # Ok(())
//! # }
//! ```
//!
//! ## 2D Interpolation
//!
//! ```
//! # use rsl_interpolation::*;
//! # fn main() -> Result<(), InterpolationError>{
//! let xa = [0.0, 1.0, 2.0, 3.0];
//! let ya = [0.0, 2.0, 4.0, 6.0];
//! // z = x + y, in column-major order
//! let za = [
//!     0.0, 1.0, 2.0, 3.0,
//!     2.0, 3.0, 4.0, 5.0,
//!     4.0, 5.0, 6.0, 7.0,
//!     6.0, 7.0, 8.0, 9.0,
//! ];
//! let interp = BicubicInterpolator::build(&xa, &ya, &za)?;
//!
//! let (x, y) = (1.2, 2.8);
//! let acc = &mut Accelerator2d::new();
//! let z = interp.eval(&xa, &ya, &za, x, y, acc)?;
//! let dzdx = interp.eval_deriv_x(&xa, &ya, &za, x, y, acc)?;
//! let dzdy = interp.eval_deriv_y(&xa, &ya, &za, x, y, acc)?;
//! let dzdxx = interp.eval_deriv_xx(&xa, &ya, &za, x, y, acc)?;
//! let dzdyy = interp.eval_deriv_yy(&xa, &ya, &za, x, y, acc)?;
//! let dzdxy = interp.eval_deriv_xy(&xa, &ya, &za, x, y, acc)?;
//! # Ok(())
//! # }
//! ```
//!
//! # Higher level interface
//!
//! + Interpolation type can be determined at runtime.
//! + The [`Spline`]/[`Spline2d`] objects own the data they was constructed with, and provide
//! the same evaluation methods as the lower-level [`Interpolation`]/[`Interpolation2d`] object,
//! without needing to provide the data arrays in every call.
//!
//! ## 1D Interpolation
//!
//! ```
//! # use rsl_interpolation::*;
//! # fn main() -> Result<(), InterpolationError>{
//! let xa = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
//! let ya = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
//! let spline = Spline::build::<CubicInterpolator>(&xa, &ya)?;
//!
//! let x = 3.9;
//! let acc = &mut Accelerator::new();
//! let y = spline.eval(x, acc)?;
//! let dy = spline.eval_deriv(x, acc)?;
//! let dy2 = spline.eval_deriv2(x, acc)?;
//! let int = spline.eval_integ(0.4, 3.1, acc)?;
//! # Ok(())
//! # }
//! ```
//!
//! ## 2D Interpolation
//!
//! ```
//! # use rsl_interpolation::*;
//! # fn main() -> Result<(), InterpolationError>{
//! let xa = [0.0, 1.0, 2.0, 3.0];
//! let ya = [0.0, 2.0, 4.0, 6.0];
//! // z = x + y, in column-major order
//! let za = [
//!     0.0, 1.0, 2.0, 3.0,
//!     2.0, 3.0, 4.0, 5.0,
//!     4.0, 5.0, 6.0, 7.0,
//!     6.0, 7.0, 8.0, 9.0,
//! ];
//! let spline = Spline2d::build::<BicubicInterpolator>(&xa, &ya, &za)?;
//!
//! let (x, y) = (1.2, 2.8);
//! let acc = &mut Accelerator2d::new();
//! let z = spline.eval(x, y, acc)?;
//! let dzdx = spline.eval_deriv_x(x, y, acc)?;
//! let dzdy = spline.eval_deriv_y(x, y, acc)?;
//! let dzdxx = spline.eval_deriv_xx(x, y, acc)?;
//! let dzdyy = spline.eval_deriv_yy(x, y, acc)?;
//! let dzdxy = spline.eval_deriv_xy(x, y, acc)?;
//! # Ok(())
//! # }
//! ```
//!
#![expect(rustdoc::broken_intra_doc_links, reason = "simple markdown")]
#![doc = include_str!("../gsl_features.md")]
//!
//! # Extra features
//!
//! + [`Accelerator2d`]: 2D Index Look-up accelerator. The generalization of the [`Accelerator`]
//! object for two dimensional interpolation.
//!
//! [`GSL's Interpolation`]: https://www.gnu.org/software/gsl/doc/html/interp.html
//! [`ndarray_linalg features`]: https://github.com/rust-ndarray/ndarray-linalg?tab=readme-ov-file#backend-features

mod accel;
mod accel2d;
mod error;
mod interp;
mod interp2d;
mod types;

mod spline;
mod spline2d;

pub(crate) use types::{
    check_if_inbounds, check_if_inbounds2d, check1d_data, check2d_data, diff, integ_eval,
};

// Public API

pub use accel::Accelerator;
pub use accel2d::Accelerator2d;

pub use error::*;
pub use interp::{BuildInterpolator, Interpolation};
pub use interp2d::{BuildInterpolator2d, Interpolation2d, z_get, z_idx, z_set};
pub use spline::Spline;
pub use spline2d::Spline2d;

pub use types::{
    AkimaInterpolator, AkimaPeriodicInterpolator, BicubicInterpolator, BilinearInterpolator,
    CubicInterpolator, CubicPeriodicInterpolator, LinearInterpolator, SteffenInterpolator,
};
