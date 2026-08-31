//! Definition of the higher level `Spline` object.

use crate::{Accelerator, Domain1dError, Interpolation};

/// 1D Higher level interface.
///
/// A Spline owns the data it is constructed with, and provides the same evaluation methods as the
/// lower-level Interpolator object, without needing to provide the data arrays in every call.
#[doc(alias = "gsl_spline")]
pub struct Spline {
    /// The lower-level [`Interpolation`] object.
    interpolator: Box<dyn Interpolation>,
    /// The owned x data.
    xa: Box<[f64]>,
    /// The owned y data.
    ya: Box<[f64]>,
}

impl Spline {
    /// Constructs a `Spline` from a [`Box<dyn Interpolation>`] and its `xa`, `ya` arrays.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut acc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    ///
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let spline = Spline::new(Box::new(interp), &xa, &ya);
    /// # Ok(())
    /// # }
    /// ```
    #[doc(alias = "gsl_spline_init")]
    #[must_use]
    pub fn new(interpolator: Box<dyn Interpolation>, xa: &[f64], ya: &[f64]) -> Self {
        Self {
            xa: xa.to_owned().into_boxed_slice(),
            ya: ya.to_owned().into_boxed_slice(),
            interpolator,
        }
    }

    /// Returns a reference to the `Spline`'s `xa` data.
    #[must_use]
    pub fn xa(&self) -> &[f64] {
        &self.xa
    }

    /// Returns a reference to the `Spline`'s `ya` data.
    #[must_use]
    pub fn ya(&self) -> &[f64] {
        &self.ya
    }

    /// Returns the interpolated value `y` for a given point `x`, using the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut acc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    ///
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let spline = Spline::new(Box::new(interp), &xa, &ya);
    ///
    /// let y = spline.eval(1.5, &mut acc)?;
    /// assert_relative_eq!(y, 3.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain1dError`] if `x` is outside the range of `xa`.
    #[doc(alias = "gsl_spline_eval")]
    #[doc(alias = "gsl_spline_eval_e")]
    pub fn eval(&self, x: f64, acc: &mut Accelerator) -> Result<f64, Domain1dError> {
        self.interpolator.eval(&self.xa, &self.ya, x, acc)
    }

    /// Returns the derivative `dy/dx` of an interpolated function for a given point `x`, using the
    /// [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut acc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    ///
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let spline = Spline::new(Box::new(interp), &xa, &ya);
    ///
    /// let dydx = spline.eval_deriv(1.5, &mut acc)?;
    /// assert_relative_eq!(dydx, 2.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain1dError`] if `x` is outside the range of `xa`.
    #[doc(alias = "gsl_spline_eval_deriv")]
    #[doc(alias = "gsl_spline_eval_deriv_e")]
    pub fn eval_deriv(&self, x: f64, acc: &mut Accelerator) -> Result<f64, Domain1dError> {
        self.interpolator.eval_deriv(&self.xa, &self.ya, x, acc)
    }

    /// Returns the second derivative `d²y/dx²` of an interpolated function for a given point `x`, using the
    /// [`Accelerator`] `acc`.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut acc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    ///
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let spline = Spline::new(Box::new(interp), &xa, &ya);
    ///
    /// let dydx = spline.eval_deriv2(1.5, &mut acc)?;
    /// assert_relative_eq!(dydx, 0.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain1dError`] if `x` is outside the range of `xa`.
    #[doc(alias = "gsl_spline_eval_deriv2")]
    #[doc(alias = "gsl_spline_eval_deriv2_e")]
    pub fn eval_deriv2(&self, x: f64, acc: &mut Accelerator) -> Result<f64, Domain1dError> {
        self.interpolator.eval_deriv2(&self.xa, &self.ya, x, acc)
    }

    /// Returns the numerical integral of an interpolated function over the range [`a` ,`b`], using the
    /// [`Accelerator`] `acc`.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut acc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    ///
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let spline = Spline::new(Box::new(interp), &xa, &ya);
    ///
    /// let int = spline.eval_integ(0.0, 2.0, &mut acc)?;
    /// assert_relative_eq!(int, 4.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain1dError`] if `a` or `b` is outside the range of xa.
    #[doc(alias = "gsl_spline_eval_integ")]
    #[doc(alias = "gsl_spline_eval_integ_e")]
    pub fn eval_integ(&self, a: f64, b: f64, acc: &mut Accelerator) -> Result<f64, Domain1dError> {
        self.interpolator.eval_integ(&self.xa, &self.ya, a, b, acc)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::*;

    #[test]
    fn arrays() {
        let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ya = [0.0, 2.0, 4.0, 6.0, 8.0];

        let interp: Box<dyn Interpolation> = Box::new(LinearInterpolator::build(&xa, &ya).unwrap());
        let spline = Spline::new(interp, &xa, &ya);

        assert_eq!(spline.xa(), xa);
        assert_eq!(spline.ya(), ya);
    }
}
