use crate::Accelerator;
use crate::DomainError;
use crate::InterpType;
use crate::Interpolation;
use crate::InterpolationError;

/// 1D Higher level interface.
///
/// A Spline owns the data it is constructed with, and provides the same evaluation methods as the
/// lower-level Interpolator object, without needing to provide the data arrays in every call.
///
/// # Example
/// ```
/// # use rsl_interpolation::*;
/// #
/// # fn main() -> Result<(), InterpolationError>{
/// let mut acc = Accelerator::new();
///
/// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
/// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
///
/// let interp = Cubic.build(&xa, &ya)?;
///
/// let typ = Cubic;
/// let spline = Spline::new(typ, &xa, &ya)?;
///
/// let x = 1.5;
/// let y_interp = interp.eval(&xa, &ya, x, &mut acc)?;
/// let y_spline = spline.eval(x, &mut acc)?;
///
/// assert_eq!(y_interp, y_spline);
/// #
/// # Ok(())
/// # }
/// ```
pub struct Spline<I, T>
where
    I: InterpType<T>,
{
    /// The lower-level [`Interpolator`].
    ///
    /// [`Interpolator`]: Interpolation#implementors
    pub interp: I::Interpolation,
    /// The owned x data.
    pub xa: Box<[T]>,
    /// The owned y data.
    pub ya: Box<[T]>,
    name: Box<str>,
    min_size: usize,
}

impl<I, T> Spline<I, T>
where
    I: InterpType<T>,
{
    /// Constructs a Spline of an Interpolation type `typ` from the data arrays `xa` and `ya`.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let typ = Cubic;
    ///
    /// let spline = Spline::new(typ, &xa, &ya)?;
    /// #
    /// # Ok(())
    /// # }
    #[doc(alias = "gsl_spline_init")]
    pub fn new(typ: I, xa: &[T], ya: &[T]) -> Result<Self, InterpolationError>
    where
        T: Clone,
    {
        Ok(Self {
            interp: typ.build(&xa, &ya)?,
            xa: xa.into(),
            ya: ya.into(),
            name: typ.name().into(),
            min_size: typ.min_size(),
        })
    }

    /// Returns the interpolated value `y` for a given point `x`, using the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut acc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let typ = Cubic;
    /// let spline = Spline::new(typ, &xa, &ya)?;
    /// #
    /// let y = spline.eval(1.5, &mut acc)?;
    ///
    /// assert_eq!(y, 3.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa`.
    #[doc(alias = "gsl_spline_eval")]
    #[doc(alias = "gsl_spline_eval_e")]
    pub fn eval(&self, x: T, acc: &mut Accelerator) -> Result<T, DomainError> {
        self.interp.eval(&self.xa, &self.ya, x, acc)
    }

    /// Returns the derivative `dy/dx` of an interpolated function for a given point `x`, using the
    /// [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut acc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let typ = Cubic;
    /// let spline = Spline::new(typ, &xa, &ya)?;
    ///
    /// let dydx = spline.eval_deriv(1.5, &mut acc)?;
    ///
    /// assert_eq!(dydx, 2.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa`.
    #[doc(alias = "gsl_spline_eval_deriv")]
    #[doc(alias = "gsl_spline_eval_deriv_e")]
    pub fn eval_deriv(&self, x: T, acc: &mut Accelerator) -> Result<T, DomainError> {
        self.interp.eval_deriv(&self.xa, &self.ya, x, acc)
    }

    /// Returns the second derivative `d²y/dx²` of an interpolated function for a given point `x`, using the
    /// [`Accelerator`] `acc`.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut acc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let typ = Cubic;
    /// let spline = Spline::new(typ, &xa, &ya)?;
    ///
    /// let dydx = spline.eval_deriv2(1.5, &mut acc)?;
    ///
    /// assert_eq!(dydx, 0.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa`.
    #[doc(alias = "gsl_spline_eval_deriv2")]
    #[doc(alias = "gsl_spline_eval_deriv2_e")]
    pub fn eval_deriv2(&self, x: T, acc: &mut Accelerator) -> Result<T, DomainError> {
        self.interp.eval_deriv2(&self.xa, &self.ya, x, acc)
    }

    #[allow(rustdoc::broken_intra_doc_links)]
    /// Returns the numerical integral of an interpolated function over the range [`a` ,`b`], using the
    /// [`Accelerator`] `acc`.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut acc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let typ = Cubic;
    /// let spline = Spline::new(typ, &xa, &ya)?;
    ///
    /// let int = spline.eval_integ(0.0, 2.0, &mut acc)?;
    ///
    /// assert_eq!(int, 4.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `a` or `b` is outside the range of xa.
    #[doc(alias = "gsl_spline_eval_integ")]
    #[doc(alias = "gsl_spline_eval_integ_e")]
    pub fn eval_integ(&self, a: T, b: T, acc: &mut Accelerator) -> Result<T, DomainError> {
        self.interp.eval_integ(&self.xa, &self.ya, a, b, acc)
    }

    /// Returns the name of the Interpolator.
    #[doc(alias = "gsl_spline_name")]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the minimum number of points required by the Interpolator.
    #[doc(alias = "gsl_spline_min_size")]
    pub fn min_size(&self) -> usize {
        self.min_size
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::Cubic;

    #[test]
    fn test_spline_creation() {
        let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ya = [0.0, 2.0, 4.0, 6.0, 8.0];

        let _spline = Spline::new(Cubic, &xa, &ya).unwrap();
    }

    #[test]
    fn test_spline_eval() {
        let xa = [0.0, 1.0, 2.0];
        let ya = [0.0, 1.0, 2.0];
        let spline = Spline::new(Cubic, &xa, &ya).unwrap();
        let mut acc = Accelerator::new();

        let x = 0.5;
        let y = spline.eval(x, &mut acc).unwrap();
        let dy = spline.eval_deriv(x, &mut acc).unwrap();
        let dy2 = spline.eval_deriv2(x, &mut acc).unwrap();
        let int = spline.eval_integ(0.0, x, &mut acc).unwrap();

        assert_eq!(y, 0.5);
        assert_eq!(dy, 1.0);
        assert_eq!(dy2, 0.0);
        assert_eq!(int, 0.125);
    }
}
