//! InterpType and Interpolation definitions

use crate::{Accelerator, DynInterpType};
use crate::{DomainError, InterpolationError};

/// Representation of an Interpolation Type.
pub trait InterpType<T> {
    /// The returned Interpolator, containing the calculated coefficients and providing the
    /// evaluation methods.
    type Interpolation: Interpolation<T> + Send + Sync;

    /// Creates an Interpolator from the data arrays `xa` and `ya`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    ///
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Cubic.build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    fn build(&self, xa: &[T], ya: &[T]) -> Result<Self::Interpolation, InterpolationError>;

    /// Returns the name of the Interpolator.
    #[doc(alias = "gsl_interp_name")]
    fn name(&self) -> &str;

    /// Returns the minimum number of points required by the Interpolator.
    #[doc(alias = "gsl_interp_min_size")]
    fn min_size(&self) -> usize;
}

/// Defines the required evaluation methods.
pub trait Interpolation<T> {
    /// Returns the interpolated value `y` for a given point `x`, using the data arrays `xa` and `ya` and
    /// the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Cubic.build(&xa, &ya)?;
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
    #[doc(alias = "gsl_interp_eval")]
    #[doc(alias = "gsl_interp_eval_e")]
    fn eval(&self, xa: &[T], ya: &[T], x: T, acc: &mut Accelerator) -> Result<T, DomainError>;

    /// Returns the derivative `dy/dx` of an interpolated function for a given point `x`, using the
    /// data arrays `xa` and `ya` and the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Cubic.build(&xa, &ya)?;
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
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Cubic.build(&xa, &ya)?;
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
    /// Returns the numerical integral of an interpolated function over the range [`a` ,`b`], using the
    /// data arrays `xa` and `ya` and the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Cubic.build(&xa, &ya)?;
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
    /// Returns a [`DomainError`] if `a` or `b` is outside the range of xa.
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

/// Creates a [`DynInterpType`] of `typ` type.
///
/// Useful when `typ` is not known at compile time.
///
/// # Example
/// ```
/// # use rsl_interpolation::*;
/// #
/// # fn main() -> Result<(), InterpolationError> {
/// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
/// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
/// let typ = "cubic";
///
/// let interp_type = make_interp_type(typ)?;
/// let interp = interp_type.build(&xa, &ya)?;
/// # Ok(())
/// # }
/// ```
pub fn make_interp_type<T>(typ: &str) -> Result<DynInterpType<T>, InterpolationError>
where
    T: crate::Num + ndarray_linalg::Lapack,
{
    use crate::*;

    match typ.to_lowercase().as_str() {
        "linear" => Ok(DynInterpType::new(Linear)),
        "cubic" => Ok(DynInterpType::new(Cubic)),
        "cubicperiodic" | "cubic periodic" => Ok(DynInterpType::new(CubicPeriodic)),
        "akima" => Ok(DynInterpType::new(Akima)),
        "akimaperiodic" | "akima periodic" => Ok(DynInterpType::new(AkimaPeriodic)),
        "steffen" => Ok(DynInterpType::new(Akima)),
        _ => Err(InterpolationError::InvalidType(typ.into())),
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::*;

    #[test]
    fn test_dyn_interp_type() {
        let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
        let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
        let mut acc = Accelerator::new();

        let x = 0.5;
        let interp_type = DynInterpType::new(Cubic);
        let interp = interp_type.build(&xa, &ya).unwrap();
        interp.eval(&xa, &ya, x, &mut acc).unwrap();
    }
}
