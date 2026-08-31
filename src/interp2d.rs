//! `Interpolation2d` trait definition.

use crate::{Accelerator2d, Domain2dError, check_if_inbounds2d};

/// Defines the required evaluation methods.
pub trait Interpolation2d {
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x + y
    /// let za = [
    ///     0.0, 1.0, 2.0,
    ///     2.0, 3.0, 4.0,
    ///     4.0, 5.0, 6.0,
    /// ];
    ///
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let z = interp.eval(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(z, 4.5);
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
        Ok(self.eval_extrap(xa, ya, za, x, y, acc))
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x + y
    /// let za = [
    ///     0.0, 1.0, 2.0,
    ///     2.0, 3.0, 4.0,
    ///     4.0, 5.0, 6.0,
    /// ];
    ///
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let z = interp.eval_extrap(&xa, &ya, &za, 3.0, 6.0, &mut acc);
    /// assert_relative_eq!(z, 9.0);
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
    ) -> f64;

    /// Returns the interpolated value `d = ∂z/∂x` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    ///
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let dzdx = interp.eval_deriv_x(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(dzdx, 3.0);
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    ///
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let dzdy = interp.eval_deriv_y(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(dzdy, 6.0);
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    ///
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let dzdx2 = interp.eval_deriv_xx(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(dzdx2, 0.0); // Linear Interpolation!
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    ///
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let dzdy2 = interp.eval_deriv_yy(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(dzdy2, 0.0); // Linear Interpolation!
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    ///
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let mut acc = Accelerator2d::new();
    ///
    /// let dzdxy = interp.eval_deriv_xy(&xa, &ya, &za, 1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(dzdxy, 0.0);
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
/// `j*len(x) + i`.
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
#[must_use]
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
#[must_use]
pub fn z_get(za: &[f64], i: usize, j: usize, xlen: usize, ylen: usize) -> f64 {
    if (i >= xlen) | (j >= ylen) {
        panic!("z-index out of range")
    };

    za[z_idx(i, j, xlen, ylen)]
}
