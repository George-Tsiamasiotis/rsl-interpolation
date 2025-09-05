use std::marker::PhantomData;

use crate::Accelerator;
use crate::InterpType;
use crate::Interpolation;
use crate::{DomainError, InterpolationError};

use crate::{check_if_inbounds, check1d_data};

const MIN_SIZE: usize = 2;

/// Linear Interpolation type.
///
/// The simplest type of interpolation.
#[doc(alias = "gsl_interp_linear")]
pub struct Linear;

impl<T> InterpType<T> for Linear
where
    T: crate::Num,
{
    type Interpolator = LinearInterp<T>;

    const MIN_SIZE: usize = MIN_SIZE;
    const NAME: &str = "Linear";

    /// Constructs a Linear Interpolator.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::InterpType;
    /// # use rsl_interpolation::Interpolation;
    /// # use rsl_interpolation::InterpolationError;
    /// # use rsl_interpolation::Linear;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Linear.build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    fn build(self, xa: &[T], ya: &[T]) -> Result<LinearInterp<T>, InterpolationError> {
        check1d_data(xa, ya, MIN_SIZE)?;
        Ok(LinearInterp {
            _variable_type: PhantomData,
        })
    }
}

// ===============================================================================================

/// Linear Interpolator.
///
/// Provides all the evaluation methods.
///
/// Should be constructed through the [`Linear`] type.
pub struct LinearInterp<T> {
    _variable_type: PhantomData<T>,
}

impl<T> Interpolation<T> for LinearInterp<T>
where
    T: crate::Num,
{
    fn eval(&self, xa: &[T], ya: &[T], x: T, acc: &mut Accelerator) -> Result<T, DomainError> {
        check_if_inbounds(xa, x)?;
        let index = acc.find(xa, x);

        let xlo = xa[index];
        let xhi = xa[index + 1];
        let ylo = ya[index];
        let yhi = ya[index + 1];

        let dx = xhi - xlo;

        debug_assert!(dx > T::zero());
        Ok(ylo + (x - xlo) / dx * (yhi - ylo))
    }

    fn eval_deriv(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        check_if_inbounds(xa, x)?;
        let index = acc.find(xa, x);

        let xlo = xa[index];
        let xhi = xa[index + 1];
        let ylo = ya[index];
        let yhi = ya[index + 1];

        let dx = xhi - xlo;
        let dy = yhi - ylo;

        debug_assert!(dx > T::zero());
        Ok(dy / dx)
    }

    /// Always returns `0`.
    #[allow(unused_variables)]
    fn eval_deriv2(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        check_if_inbounds(xa, x)?;
        Ok(T::zero())
    }

    fn eval_integ(
        &self,
        xa: &[T],
        ya: &[T],
        a: T,
        b: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        check_if_inbounds(xa, a)?;
        check_if_inbounds(xa, b)?;
        let index_a = acc.find(xa, a);
        let index_b = acc.find(xa, b);

        let mut result = T::zero();

        for i in index_a..=index_b {
            let xlo = xa[i];
            let xhi = xa[i + 1];
            let ylo = ya[i];
            let yhi = ya[i + 1];

            let dx = xhi - xlo;
            let d = (yhi - ylo) / dx;
            let half = T::from(0.5).unwrap();

            if dx.is_zero() {
                continue;
            }

            if (i == index_a) | (i == index_b) {
                let x1 = if i == index_a { a } else { xlo };
                let x2 = if i == index_b { b } else { xhi };
                result += (x2 - x1) * (ylo + half * d * ((x2 - xlo) + (x1 - xlo)));
            } else {
                result += half * dx * (ylo + yhi);
            }
        }
        Ok(result)
    }
}
