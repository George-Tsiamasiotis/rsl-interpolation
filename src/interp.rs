//! `Interpolation` and `BuildInterpolator` traits definition.

use crate::{Accelerator, Domain1dError, InterpolationError};

/// 1D Interpolator build method.
pub trait BuildInterpolator: Interpolation + Sized {
    /// The minimum required number of data points.
    #[doc(alias = "gsl_interp_min_size")]
    const MIN_SIZE: usize;

    /// Builds the Interpolator.
    #[doc(alias = "gsl_interp_init")]
    #[expect(clippy::missing_errors_doc, reason = "documented on the implementors")]
    fn build(xa: &[f64], ya: &[f64]) -> Result<Self, InterpolationError>;
}

/// Defines the available interpolation methods.
pub trait Interpolation: Send + Sync + 'static {
    /// Returns the interpolated value `y` for a given point `x`, using the data arrays `xa` and `ya` and
    /// the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let mut acc = Accelerator::new();
    ///
    /// let y = interp.eval(&xa, &ya, 1.5, &mut acc)?;
    /// assert_relative_eq!(y, 3.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain1dError`] if `x` is outside the range of `xa`.
    #[doc(alias = "gsl_interp_eval")]
    #[doc(alias = "gsl_interp_eval_e")]
    fn eval(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError>;

    /// Returns the derivative `dy/dx` of an interpolated function for a given point `x`, using the
    /// data arrays `xa` and `ya` and the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let mut acc = Accelerator::new();
    ///
    /// let dydx = interp.eval_deriv(&xa, &ya, 1.5, &mut acc)?;
    /// assert_relative_eq!(dydx, 2.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain1dError`] if `x` is outside the range of `xa`.
    #[doc(alias = "gsl_interp_eval_deriv")]
    #[doc(alias = "gsl_interp_eval_deriv_e")]
    fn eval_deriv(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError>;

    /// Returns the second derivative `d²y/dx²` of an interpolated function for a given point `x`, using the
    /// data arrays `xa` and `ya` and the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let mut acc = Accelerator::new();
    ///
    /// let dydx2 = interp.eval_deriv2(&xa, &ya, 1.5, &mut acc)?;
    /// assert_relative_eq!(dydx2, 0.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain1dError`] if `x` is outside the range of `xa`.
    #[doc(alias = "gsl_interp_eval_deriv2")]
    #[doc(alias = "gsl_interp_eval_deriv2_e")]
    fn eval_deriv2(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError>;

    /// Returns the numerical integral of an interpolated function over the range [`a` ,`b`], using the
    /// data arrays `xa` and `ya` and the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let mut acc = Accelerator::new();
    ///
    /// let int = interp.eval_integ(&xa, &ya, 0.0, 2.0, &mut acc)?;
    /// assert_relative_eq!(int, 4.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain1dError`] if `a` or `b` is outside the range of xa.
    #[doc(alias = "gsl_interp_eval_integ")]
    #[doc(alias = "gsl_interp_eval_integ_e")]
    fn eval_integ(
        &self,
        xa: &[f64],
        ya: &[f64],
        a: f64,
        b: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError>;
}
