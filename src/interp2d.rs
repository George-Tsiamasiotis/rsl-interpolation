//! `Interpolator2d`, `Interpolation2d` and `DynInterpolator2d` trait definitions.
//!
//! NOTE: `Interpolator2d` must be a separate trait for dyn compatibility.

use std::fmt::Debug;

use crate::types::check_if_inbounds2d;
use crate::{Accelerator2d, Domain2dError, InterpolatorError};
use crate::{BicubicInterpolator, BilinearInterpolator};

/// Available 2D Interpolation Types.
///
/// ## References
///
/// Numerical Algorithms with C - Gisela Engeln-Mullges, Frank Uhlig - 1996 -
/// Algorithm 10.1, pg 254
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interpolation2dType {
    /// Bilinear interpolation.
    ///
    /// This interpolation method does not require any additional memory.
    Bilinear,
    /// Bicubic Interpolation.
    Bicubic,
}

/// Defines the Interpolator build method.
///
/// > # **Important**
/// >
/// > The `za` array must be defined in **column-major (Fortran)** style. This is done to comply
/// > with GSL's interface.
/// >
pub trait Interpolator2d: Clone + Sized {
    /// Creates a 2D Interpolator from the data arrays `xa`, `ya` and `za`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0, 3.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0];
    /// // z = x + y
    /// let za = [
    ///     0.0, 1.0, 2.0, 3.0,
    ///     2.0, 3.0, 4.0, 5.0,
    ///     4.0, 5.0, 6.0, 7.0,
    ///     6.0, 7.0, 8.0, 9.0,
    /// ];
    ///
    /// let interp = BicubicInterpolator::build(&xa, &ya, &za)?;
    /// # Ok(())
    /// # }
    /// ```
    #[doc(alias = "gsl_interp2d_init")]
    fn build(xa: &[f64], ya: &[f64], za: &[f64]) -> Result<Self, InterpolatorError>;

    /// Returns the name of the Interpolator.
    #[doc(alias = "gsl_interp2d_name")]
    fn name(&self) -> &str;

    /// Returns the minimum number of points required by the Interpolator.
    #[doc(alias = "gsl_interp2d_min_size")]
    fn min_size(&self) -> usize;
}

/// Defines the required evaluation methods.
#[allow(private_bounds)]
pub trait Interpolation2d: DynInterpolator2dClone + Send + Sync + 'static + Debug {
    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerator2d`] `acc`.
    ///
    /// # Note
    ///
    /// This function only performs the bounds check, and then calls `eval_extrap()`, where the
    /// actual evaluation is implemented.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x + y
    /// let za = [
    ///     0.0, 1.0, 2.0,
    ///     2.0, 3.0, 4.0,
    ///     4.0, 5.0, 6.0,
    /// ];
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let z = interp.eval(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(z, 4.5);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside
    /// the range of `ya`.
    #[doc(alias = "gsl_interp2d_eval")]
    #[doc(alias = "gsl_interp2d_eval_e")]
    fn eval(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        // Calculation is the same, with the added bounds check
        check_if_inbounds2d(xa, ya, x, y)?;
        self.eval_extrap(xa, ya, za, x, y, acc)
    }

    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerator2d`] `acc`.
    ///
    /// # Note
    ///
    /// This function performs *no bound checking*, so when `x` is outside the range of `xa` or y
    /// is outside the range of `ya`, extrapolation is performed.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x + y
    /// let za = [
    ///     0.0, 1.0, 2.0,
    ///     2.0, 3.0, 4.0,
    ///     4.0, 5.0, 6.0,
    /// ];
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let z = interp.eval_extrap(&xa, &ya, &za, 3.0, 6.0, &mut acc)?;
    ///
    /// assert_eq!(z, 9.0);
    /// # Ok(())
    /// # }
    /// ```
    #[doc(alias = "gsl_interp2d_eval_extrap")]
    #[doc(alias = "gsl_interp2d_eval_extrap_e")]
    fn eval_extrap(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError>;

    /// Returns the interpolated value `d = ∂z/∂x` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let dzdx = interp.eval_deriv_x(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(dzdx, 3.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_interp2d_eval_deriv_x")]
    #[doc(alias = "gsl_interp2d_eval_deriv_x_e")]
    fn eval_deriv_x(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError>;

    /// Returns the interpolated value `d = ∂z/∂y` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let dzdy = interp.eval_deriv_y(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(dzdy, 6.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_interp2d_eval_deriv_y")]
    #[doc(alias = "gsl_interp2d_eval_deriv_y_e")]
    fn eval_deriv_y(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError>;

    /// Returns the interpolated value `d = 𝜕²z/𝜕x²` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let dzdx2 = interp.eval_deriv_xx(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(dzdx2, 0.0); // Linear Interpolation!
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_interp2d_eval_deriv_xx")]
    #[doc(alias = "gsl_interp2d_eval_deriv_xx_e")]
    fn eval_deriv_xx(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError>;

    /// Returns the interpolated value `d = 𝜕²z/𝜕y²` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let dzdy2 = interp.eval_deriv_yy(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(dzdy2, 0.0); // Linear Interpolation!
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_interp2d_eval_deriv_yy")]
    #[doc(alias = "gsl_interp2d_eval_deriv_yy_e")]
    fn eval_deriv_yy(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError>;

    /// Returns the interpolated value `d = 𝜕²z/𝜕x𝜕y` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let dzdxy = interp.eval_deriv_xy(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(dzdxy, 0.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_interp2d_eval_deriv_xy")]
    #[doc(alias = "gsl_interp2d_eval_deriv_xy_e")]
    fn eval_deriv_xy(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError>;
}

/// Returns the index corresponding to the grid point (`i`, `j`). The index is given by
/// `j*len(x) + i`
///
/// # Important
///
/// The `za` array is indexed in column-major style (Fortran style), so it must be defined
/// accordingly.
///
/// # Example
///
/// ```
/// # use rsl_interpolation::*;
/// let xa = [0.0, 1.0];
/// let ya = [0.0, 2.0];
/// let za = [
///     0.0, 1.0, // <- This one
///     2.0, 3.0,
/// ];
/// let za_index = z_idx(0, 1, xa.len(), ya.len());
/// assert_eq!(za_index, 2);
/// ```
///
/// # Panics
///
/// Panics if `i>=xlen` or `i>=ylen`.
#[doc(alias = "gsl_interp2d_idx")]
pub fn z_idx(i: usize, j: usize, xlen: usize, ylen: usize) -> usize {
    if (i >= xlen) | (j >= ylen) {
        panic!("z-index out of range")
    } else {
        j * xlen + i
    }
}

/// Sets the value `z` of grid point (`i`, `j`) of the array `za` to `z`.
///
/// # Important
///
/// The `za` array is indexed in column-major style (Fortran style), so it must be defined
/// accordingly.
///
/// # Example
///
/// ```
/// # use rsl_interpolation::*;
/// let xa = [0.0, 1.0];
/// let ya = [0.0, 2.0];
/// let mut za = [
///     0.0, 1.0, // <- We set this one
///     2.0, 3.0,
/// ];
/// z_set(&mut za, 10.0, 0, 1, xa.len(), ya.len());
/// assert_eq!(za[2], 10.0);
/// ```
///
/// # Panics
///
/// Panics if `i>=xlen` or `j>=ylen`.
#[doc(alias = "gsl_inter2d_set")]
pub fn z_set<T>(za: &mut [T], z: T, i: usize, j: usize, xlen: usize, ylen: usize) {
    if (i >= xlen) | (j >= ylen) {
        panic!("z-index out of range")
    };

    za[z_idx(i, j, xlen, ylen)] = z;
}

/// Returns the value `z` of grid point (`i`, `j`) of the array `za` to `z`.
///
/// # Important
///
/// The `za` array is indexed in column-major style (Fortran style), so it must be defined
/// accordingly.
///
/// # Example
///
/// ```
/// # use rsl_interpolation::*;
/// let xa = [0.0, 1.0];
/// let ya = [0.0, 2.0];
/// let za = [
///     0.0, 10.0, // <- We want this one
///     2.0, 3.0,
/// ];
/// let g = z_get(&za, 1, 0, xa.len(), ya.len());
/// assert_eq!(g, 10.0);
/// ```
///
/// # Panics
///
/// Panics if `i>=xlen` or `j>=ylen`.
#[doc(alias = "gsl_inter2d_get")]
pub fn z_get(za: &[f64], i: usize, j: usize, xlen: usize, ylen: usize) -> f64 {
    if (i >= xlen) | (j >= ylen) {
        panic!("z-index out of range")
    };

    za[z_idx(i, j, xlen, ylen)]
}

/// Interpolator with dynamically determined type.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct DynInterpolator2d {
    interpolator: Box<dyn Interpolation2d>,
    typ: Interpolation2dType,
}

impl DynInterpolator2d {
    /// Creates a `DynInterpolator` of `typ` type.
    ///
    /// Useful when `typ` is not known at compile time.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError> {
    /// let xa = [0.0, 1.0, 2.0, 3.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0];
    /// // z = x + y, in column-major order
    /// let za = [
    ///     0.0, 1.0, 2.0, 3.0,
    ///     2.0, 3.0, 4.0, 5.0,
    ///     4.0, 5.0, 6.0, 7.0,
    ///     6.0, 7.0, 8.0, 9.0,
    /// ];
    /// let typ2d = Interpolation2dType::Bicubic;
    ///
    /// let interp2d = DynInterpolator2d::build(typ2d, &xa, &ya, &za)?;
    ///
    /// assert_eq!(interp2d.typ(), typ2d);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`InterpolatorError`] if the Interpolator fails to build.
    pub fn build(
        typ: Interpolation2dType,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
    ) -> Result<Self, InterpolatorError> {
        use Interpolation2dType::*;
        let interpolator: Box<dyn Interpolation2d> = match typ {
            Bilinear => Box::new(BilinearInterpolator::build(xa, ya, za)?),
            Bicubic => Box::new(BicubicInterpolator::build(xa, ya, za)?),
        };
        Ok(Self { interpolator, typ })
    }

    /// Returns the interpolator's [`Interpolation2dType`].
    pub fn typ(&self) -> Interpolation2dType {
        self.typ
    }
}

impl DynInterpolator2d {
    /// Creates a `DynInterpolator2d` of `typ` type.
    ///
    /// Useful when `typ` is not known at compile time.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError> {
    /// let xa = [0.0, 1.0, 2.0, 3.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0];
    /// // z = x + y, in column-major order
    /// let za = [
    ///     0.0, 1.0, 2.0, 3.0,
    ///     2.0, 3.0, 4.0, 5.0,
    ///     4.0, 5.0, 6.0, 7.0,
    ///     6.0, 7.0, 8.0, 9.0,
    /// ];
    /// let typ2d = Interpolation2dType::Bicubic;
    ///
    /// let interp2d = DynInterpolator2d::build(typ2d, &xa, &ya, &za)?;
    ///
    /// assert_eq!(interp2d.typ(), typ2d);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`InterpolatorError`] if the Interpolator fails to build.
    pub fn new(
        typ: Interpolation2dType,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
    ) -> Result<Self, InterpolatorError> {
        use Interpolation2dType::*;
        let interpolator: Box<dyn Interpolation2d> = match typ {
            Bilinear => Box::new(BilinearInterpolator::build(xa, ya, za)?),
            Bicubic => Box::new(BicubicInterpolator::build(xa, ya, za)?),
        };
        Ok(Self { interpolator, typ })
    }
}

impl Interpolation2d for DynInterpolator2d {
    fn eval_extrap(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator.eval_extrap(xa, ya, za, x, y, acc)
    }

    fn eval_deriv_x(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator.eval_deriv_x(xa, ya, za, x, y, acc)
    }

    fn eval_deriv_y(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator.eval_deriv_y(xa, ya, za, x, y, acc)
    }

    fn eval_deriv_xx(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator.eval_deriv_xx(xa, ya, za, x, y, acc)
    }

    fn eval_deriv_yy(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator.eval_deriv_yy(xa, ya, za, x, y, acc)
    }

    fn eval_deriv_xy(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator.eval_deriv_xy(xa, ya, za, x, y, acc)
    }
}

// HACK: make `DynInterpolator` Clone
// https://stackoverflow.com/questions/30353462/how-to-clone-a-struct-storing-a-boxed-trait-object

trait DynInterpolator2dClone {
    fn clone_box(&self) -> Box<dyn Interpolation2d>;
}

impl<T> DynInterpolator2dClone for T
where
    T: 'static + Interpolation2d + Clone,
{
    fn clone_box(&self) -> Box<dyn Interpolation2d> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Interpolation2d> {
    fn clone(&self) -> Box<dyn Interpolation2d> {
        self.clone_box()
    }
}
