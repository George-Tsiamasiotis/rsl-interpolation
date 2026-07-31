//! A re-write of [`GSL's Interpolation`] in Rust.
//!
//! [`GSL's Interpolation`]: https://www.gnu.org/software/gsl/doc/html/interp.html
//!
//! # Example - 1D lower level interface
//!
//! + Interpolation type must be known at compilation time.
//! + The [`Interpolator`] object does not store the arrays it was constructed with; they must be
//! passed as arguments at each evaluation.
//!
//! ```
//! # use rsl_interpolation::*;
//! #
//! # fn main() -> Result<(), InterpolatorError>{
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
//! # Example - 2D lower level interface
//!
//! + Interpolation type must be known at compilation time.
//! + The [`Interpolator2d`] object does not store the arrays it was constructed with; they must be
//! passed as arguments at each evaluation.
//!
//! ```
//! # use rsl_interpolation::*;
//! #
//! # fn main() -> Result<(), InterpolatorError>{
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
//! # Example - Dynamically determined interpolation type
//!
//! + Interpolation type can be determined at runtime.
//! + The [`DynInterpolator`] and [`DynInterpolator2d`] objects do not store the arrays they were
//! constructed with; they must be passed as arguments at each evaluation.
//!
//! ```
//! # use rsl_interpolation::*;
//! #
//! # fn main() -> Result<(), InterpolatorError>{
//! let xa = [0.0, 1.0, 2.0, 3.0];
//! let ya = [0.0, 2.0, 4.0, 6.0];
//! // z = x + y, in column-major order
//! let za = [
//!     0.0, 1.0, 2.0, 3.0,
//!     2.0, 3.0, 4.0, 5.0,
//!     4.0, 5.0, 6.0, 7.0,
//!     6.0, 7.0, 8.0, 9.0,
//! ];
//!
//! // 1D Interpolation
//! let typ = Interpolation1dType::Cubic;
//! let dyn_interp = DynInterpolator::build(typ, &xa, &ya)?;
//!
//! let x = 1.8;
//! let acc = &mut Accelerator::new();
//! let y = dyn_interp.eval(&xa, &ya, x, acc)?;
//!
//! // 2D Interpolation
//! let typ2d = Interpolation2dType::Bicubic;
//! let dyn_interp2d = DynInterpolator2d::build(typ2d, &xa, &ya, &za)?;
//!
//! let (x, y) = (1.2, 2.8);
//! let acc = &mut Accelerator2d::new();
//! let z = dyn_interp2d.eval(&xa, &ya, &za, x, y, acc)?;
//! # Ok(())
//! # }
//! ```
//!
//! # Example - 1D higher level interface
//!
//! + Interpolation type can be determined at runtime.
//! + The [`Spline`] object owns the data it was constructed with, and provides the same evaluation
//! methods as the lower-level [`Interpolation`] object, without needing to provide the data arrays
//! in every call.
//!
//! ```
//! # use rsl_interpolation::*;
//! #
//! # fn main() -> Result<(), InterpolatorError>{
//! let xa = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
//! let ya = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
//! let interp = DynInterpolator::build(Interpolation1dType::Cubic, &xa, &ya)?;
//! let spline = Spline::new(interp, &xa, &ya);
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
//! # Example - 2D higher level interface
//!
//! + Interpolation type can be determined at runtime.
//! + The [`Spline2d`] object owns the data it was constructed with, and provides the same evaluation
//! methods as the lower-level [`Interpolation2d`] object, without needing to provide the data arrays
//! in every call.
//!
//! ```
//! # use rsl_interpolation::*;
//! #
//! # fn main() -> Result<(), InterpolatorError>{
//! let xa = [0.0, 1.0, 2.0, 3.0];
//! let ya = [0.0, 2.0, 4.0, 6.0];
//! // z = x + y, in column-major order
//! let za = [
//!     0.0, 1.0, 2.0, 3.0,
//!     2.0, 3.0, 4.0, 5.0,
//!     4.0, 5.0, 6.0, 7.0,
//!     6.0, 7.0, 8.0, 9.0,
//! ];
//! let interp = DynInterpolator2d::build(Interpolation2dType::Bicubic, &xa, &ya, &za)?;
//! let spline = Spline2d::new(interp, &xa, &ya, &za);
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
//! # 1D Interpolation Types
//!
//! + [LinearInterpolator]
//! + [CubicInterpolator]
//! + [CubicPeriodicInterpolator]
//! + [AkimaInterpolator]
//! + [AkimaPeriodicInterpolator]
//! + [SteffenInterpolator]
//! + [DynInterpolator] - Dynamically determined interpolation type
//!
//! # 2D Interpolation Types
//!
//! + [BilinearInterpolator]
//! + [BicubicInterpolator]
//!
//! # Higher level Interface
//!
//! + [Spline] - Dynamically determined 1D interpolator that owns the data it was constructed with.
//! + [Spline2d] - Dynamically determined 2D interpolator that owns the data it was constructed with.
//!
#![allow(rustdoc::broken_intra_doc_links)]
#![doc = include_str!("../TODO.md")]
//!
//! # Extra features
//!
//! + [`Accelerator2d`]: 2D Index Look-up accelerator. The generalization of the [`Accelerator`]
//! object for two dimensional interpolation.

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
pub use interp::{DynInterpolator, Interpolation, Interpolation1dType, Interpolator};
pub use interp2d::{DynInterpolator2d, Interpolation2d, Interpolation2dType, Interpolator2d};
pub use interp2d::{z_get, z_idx, z_set};
pub use spline::Spline;
pub use spline2d::Spline2d;

pub use types::{
    AkimaInterpolator, AkimaPeriodicInterpolator, BicubicInterpolator, BilinearInterpolator,
    CubicInterpolator, CubicPeriodicInterpolator, LinearInterpolator, SteffenInterpolator,
};
