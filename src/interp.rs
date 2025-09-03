use crate::Accelerator;
use crate::{DomainError, InterpolationError};

/// Defines the required methods for every Interpolation type.
pub trait Interpolation<T>
where
    T: num::Float + std::fmt::Debug,
{
    /// The minimum number of points required by the interpolator. For example, Akima spline
    /// interpolation requires a minimum of 5 points.
    const MIN_SIZE: usize;

    /// The name of the interpolator.
    const NAME: &'static str;

    /// Creates a new Interpolator for the data (`xa`, `ya`), where `xa` and `ya` are slices of the
    /// x and y data points.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::Interpolation;
    /// # use rsl_interpolation::InterpolationError;
    /// # use rsl_interpolation::Cubic;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Cubic::new(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    #[doc(alias = "gsl_interp_init")]
    fn new(xa: &[T], ya: &[T]) -> Result<Self, InterpolationError>
    where
        Self: Sized;

    /// Returns the interpolated value `y` for a given point `x`, using the data arrays `xa` and `ya` and
    /// the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::Interpolation;
    /// # use rsl_interpolation::InterpolationError;
    /// # use rsl_interpolation::Accelerator;
    /// # use rsl_interpolation::Cubic;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Cubic::new(&xa, &ya)?;
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
    /// Returns a [`DomainError`] if `x` is outside the range of `xa`.
    ///
    /// [`Accelerator`]: struct.Accelerator.html
    /// [`DomainError`]: struct.DomainError.html
    #[doc(alias = "gsl_interp_eval")]
    #[doc(alias = "gsl_interp_eval_e")]
    fn eval(&self, xa: &[T], ya: &[T], x: T, acc: &mut Accelerator) -> Result<T, DomainError>;

    /// Returns the derivative `dy/dx` of an interpolated function for a given point `x`, using the
    /// data arrays `xa` and `ya` and the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::Interpolation;
    /// # use rsl_interpolation::InterpolationError;
    /// # use rsl_interpolation::Accelerator;
    /// # use rsl_interpolation::Cubic;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Cubic::new(&xa, &ya)?;
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
    /// Returns a [`DomainError`] if `x` is outside the range of `xa`.
    ///
    /// [`DomainError`]: struct.DomainError.html
    /// [`Accelerator`]: struct.Accelerator.html
    #[doc(alias = "gsl_interp_eval_deriv")]
    #[doc(alias = "gsl_interp_eval_deriv_e")]
    fn eval_deriv(&self, xa: &[T], ya: &[T], x: T, acc: &mut Accelerator)
    -> Result<T, DomainError>;

    /// Returns the second derivative `d²y/dx²` of an interpolated function for a given point `x`, using the
    /// data arrays `xa` and `ya` and the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::Interpolation;
    /// # use rsl_interpolation::InterpolationError;
    /// # use rsl_interpolation::Accelerator;
    /// # use rsl_interpolation::Cubic;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Cubic::new(&xa, &ya)?;
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
    /// Returns a [`DomainError`] if `x` is outside the range of `xa`.
    ///
    /// [`DomainError`]: struct.DomainError.html
    /// [`Accelerator`]: struct.Accelerator.html
    #[doc(alias = "gsl_interp_eval_deriv2")]
    #[doc(alias = "gsl_interp_eval_deriv2_e")]
    fn eval_deriv2(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError>;

    #[allow(rustdoc::broken_intra_doc_links)]
    /// Returns the numerical integral of an interpolated function over the range [`a`,`b`], using the
    /// data arrays `xa` and `ya` and the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::Interpolation;
    /// # use rsl_interpolation::InterpolationError;
    /// # use rsl_interpolation::Accelerator;
    /// # use rsl_interpolation::Cubic;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Cubic::new(&xa, &ya)?;
    /// let mut acc = Accelerator::new();
    ///
    /// let dydx2 = interp.eval_integ(&xa, &ya, 0.0, 2.0, &mut acc)?;
    ///
    /// assert_eq!(dydx2, 4.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `a` or `b` is outside the range of xa.
    ///
    /// [`DomainError`]: struct.DomainError.html
    /// [`Accelerator`]: struct.Accelerator.html
    #[doc(alias = "gsl_interp_eval_integ")]
    #[doc(alias = "gsl_interp_eval_integ_e")]
    fn eval_integ(
        &self,
        xa: &[T],
        ya: &[T],
        a: T,
        b: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError>;
}
