//! `Interpolator`, `Interpolation` and `DynInterpolator` trait definitions.
//!
//! NOTE: `Interpolator` must be a separate trait for dyn compatibility.

use std::fmt::Debug;

use crate::{Accelerator, Domain1dError, InterpolatorError};
use crate::{
    AkimaInterpolator, AkimaPeriodicInterpolator, CubicInterpolator, CubicPeriodicInterpolator,
    LinearInterpolator, SteffenInterpolator,
};

/// Available 1D Interpolation Types.
///
/// ## References
///
/// Numerical Algorithms with C - Gisela Engeln-Mullges, Frank Uhlig - 1996 -
/// Algorithm 10.1, pg 254
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interpolation1dType {
    /// Linear interpolation.
    ///
    /// This interpolation method does not require any additional memory.
    Linear,
    /// Cubic spline with natural boundary conditions.
    ///
    /// The resulting curve is piecewise cubic on each interval, with matching first and second
    /// derivatives at the supplied data-points. The second derivative is chosen to be zero at the
    /// first point and last point.
    Cubic,
    /// Cubic spline with periodic boundary conditions.
    ///
    /// The resulting curve is piecewise cubic on each interval, with matching first and second
    /// derivatives at the supplied data-points. The derivatives at the first and last points are
    /// also matched. Note that the last point in the data must have the same y-value as the first
    /// point, otherwise the resulting periodic interpolation will have a discontinuity at the
    /// boundary.
    CubicPeriodic,
    /// Non-rounded Akima spline with natural boundary conditions.
    ///
    /// This method uses the non-rounded corner algorithm of Wodicka.
    Akima,
    /// Non-rounded Akima spline with periodic boundary conditions.
    ///
    /// This method uses the non-rounded corner algorithm of Wodicka.
    AkimaPeriodic,
    /// Steffen interpolation.
    ///
    /// Steffen’s method guarantees the monotonicity of the interpolating function between the
    /// given data points. Therefore, minima and maxima can only occur exactly at the data
    /// points, and there can never be spurious oscillations between data points. The interpolated
    /// function is piecewise cubic in each interval. The resulting curve and its first derivative
    /// are guaranteed to be continuous, but the second derivative may be discontinuous.
    Steffen,
}

/// Defines the Interpolator build method.
pub trait Interpolator: Clone + Sized {
    /// Creates an Interpolator from the data arrays `xa` and `ya`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`InterpolatorError`] if the Interpolator fails to build.
    #[doc(alias = "gsl_interp_init")]
    fn build(xa: &[f64], ya: &[f64]) -> Result<Self, InterpolatorError>;

    /// Returns the name of the Interpolator.
    #[doc(alias = "gsl_interp_name")]
    fn name(&self) -> &str;

    /// Returns the minimum number of points required by the Interpolator.
    #[doc(alias = "gsl_interp_min_size")]
    fn min_size(&self) -> usize;
}

/// Defines the available interpolation methods.
#[allow(private_bounds)]
pub trait Interpolation: DynInterpolatorClone + Send + Sync + 'static + Debug {
    /// Returns the interpolated value `y` for a given point `x`, using the data arrays `xa` and `ya` and
    /// the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let mut acc = Accelerator::new();
    ///
    /// let y = interp.eval(&xa, &ya, 1.5, &mut acc)?;
    ///
    /// assert_eq!(y, 3.0);
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
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let mut acc = Accelerator::new();
    ///
    /// let dydx = interp.eval_deriv(&xa, &ya, 1.5, &mut acc)?;
    ///
    /// assert_eq!(dydx, 2.0);
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
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let mut acc = Accelerator::new();
    ///
    /// let dydx2 = interp.eval_deriv2(&xa, &ya, 1.5, &mut acc)?;
    ///
    /// assert_eq!(dydx2, 0.0);
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

    #[allow(rustdoc::broken_intra_doc_links)]
    /// Returns the numerical integral of an interpolated function over the range [`a` ,`b`], using the
    /// data arrays `xa` and `ya` and the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// let mut acc = Accelerator::new();
    ///
    /// let int = interp.eval_integ(&xa, &ya, 0.0, 2.0, &mut acc)?;
    ///
    /// assert_eq!(int, 4.0);
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

/// Interpolator with dynamically determined type.
#[derive(Clone, Debug)]
#[non_exhaustive]
pub struct DynInterpolator {
    interpolator: Box<dyn Interpolation>,
    typ: Interpolation1dType,
}

impl DynInterpolator {
    /// Creates a `DynInterpolator` of `typ` type.
    ///
    /// Useful when `typ` is not known at compile time.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError> {
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let typ = Interpolation1dType::Akima;
    ///
    /// let interp = DynInterpolator::build(typ, &xa, &ya)?;
    ///
    /// assert_eq!(interp.typ(), typ);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`InterpolatorError`] if the Interpolator fails to build.
    pub fn build(
        typ: Interpolation1dType,
        xa: &[f64],
        ya: &[f64],
    ) -> Result<Self, InterpolatorError> {
        use Interpolation1dType::*;
        let interpolator: Box<dyn Interpolation> = match typ {
            Linear => Box::new(LinearInterpolator::build(xa, ya)?),
            Cubic => Box::new(CubicInterpolator::build(xa, ya)?),
            CubicPeriodic => Box::new(CubicPeriodicInterpolator::build(xa, ya)?),
            Akima => Box::new(AkimaInterpolator::build(xa, ya)?),
            AkimaPeriodic => Box::new(AkimaPeriodicInterpolator::build(xa, ya)?),
            Steffen => Box::new(SteffenInterpolator::build(xa, ya)?),
        };
        Ok(Self { interpolator, typ })
    }

    /// Returns the interpolator's [`Interpolation1dType`].
    pub fn typ(&self) -> Interpolation1dType {
        self.typ
    }
}

impl Interpolation for DynInterpolator {
    fn eval(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        self.interpolator.eval(xa, ya, x, acc)
    }

    fn eval_deriv(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        self.interpolator.eval_deriv(xa, ya, x, acc)
    }

    fn eval_deriv2(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        self.interpolator.eval_deriv2(xa, ya, x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[f64],
        ya: &[f64],
        a: f64,
        b: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        self.interpolator.eval_integ(xa, ya, a, b, acc)
    }
}

// HACK: make `DynInterpolator` Clone
// https://stackoverflow.com/questions/30353462/how-to-clone-a-struct-storing-a-boxed-trait-object

trait DynInterpolatorClone {
    fn clone_box(&self) -> Box<dyn Interpolation>;
}

impl<T> DynInterpolatorClone for T
where
    T: 'static + Interpolation + Clone,
{
    fn clone_box(&self) -> Box<dyn Interpolation> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Interpolation> {
    fn clone(&self) -> Box<dyn Interpolation> {
        self.clone_box()
    }
}
