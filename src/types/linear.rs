//! Implementor for Linear Interpolator.

use std::fmt::Debug;
use std::marker::PhantomData;

use crate::Accelerator;
use crate::Interpolation;
use crate::{DomainError, InterpolationError};

use crate::types::utils::check_data;
use crate::types::utils::check_if_inbounds;

/// Linear interpolation.
///
/// The simplest type of interpolation.
///
/// ## Example
///
/// ```
/// # use rsl_interpolation::Interpolation;
/// # use rsl_interpolation::InterpolationError;
/// # use rsl_interpolation::Linear;
/// # use rsl_interpolation::Accelerator;
/// #
/// # fn main() -> Result<(), InterpolationError>{
/// let xa = [0.0, 1.0, 2.0];
/// let ya = [0.0, 2.0, 4.0];
/// let interp = Linear::new(&xa, &ya)?;
/// # Ok(())
/// # }
/// ```
pub struct Linear<T> {
    _variable_type: PhantomData<T>,
}

impl<T> Interpolation<T> for Linear<T>
where
    T: num::Float + Debug + std::ops::AddAssign,
{
    const MIN_SIZE: usize = 2;
    const NAME: &'static str = "linear";

    fn new(xa: &[T], ya: &[T]) -> Result<Self, InterpolationError> {
        check_data(xa, ya, Self::MIN_SIZE)?;
        Ok(Self {
            _variable_type: PhantomData,
        })
    }

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
