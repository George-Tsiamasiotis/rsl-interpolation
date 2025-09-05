use ndarray_linalg::Lapack;

use crate::Accelerator;
use crate::DomainError;
use crate::InterpType;
use crate::Interpolation;
use crate::InterpolationError;

#[allow(dead_code)]
pub struct Spline<I, T> {
    interp: I,
    xa: Vec<T>,
    ya: Vec<T>,
}

impl<I, T> Spline<I, T>
where
    I: Interpolation<T>,
    T: crate::Num + Lapack,
{
    /// Constructs a Spline of an Interpolation type `typ` from the data arrays `xa` and `ya`.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::Spline;
    /// # use rsl_interpolation::Cubic;
    /// # use rsl_interpolation::InterpType;
    /// # use rsl_interpolation::InterpolationError;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let typ = Cubic;
    ///
    /// let spline = Spline::build(typ, &xa, &ya)?;
    /// #
    /// # Ok(())
    /// # }
    pub fn build(
        typ: impl InterpType<T, Interpolator = I>,
        xa: &[T],
        ya: &[T],
    ) -> Result<Self, InterpolationError> {
        let xa = xa.to_owned();
        let ya = ya.to_owned();

        let interp = typ.build(&xa, &ya)?;

        let spline = Self { interp, xa, ya };

        Ok(spline)
    }

    /// Returns the interpolated value `y` for a given point `x`, using the [`Accelerator`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::Spline;
    /// # use rsl_interpolation::Cubic;
    /// # use rsl_interpolation::InterpType;
    /// # use rsl_interpolation::Accelerator;
    /// # use rsl_interpolation::InterpolationError;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut acc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let typ = Cubic;
    /// let spline = Spline::build(typ, &xa, &ya)?;
    /// #
    /// let y = spline.eval(1.5, &mut acc)?;
    ///
    /// assert_eq!(y, 3.0);
    /// # Ok(())
    /// # }
    /// ```
    pub fn eval(&self, x: T, acc: &mut Accelerator) -> Result<T, DomainError> {
        self.interp.eval(&self.xa, &self.ya, x, acc)
    }
}
