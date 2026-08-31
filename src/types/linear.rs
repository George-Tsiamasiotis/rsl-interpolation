//! Definition of `Linear` Interpolator.

use crate::{Accelerator, Domain1dError, Interpolation, InterpolationError};
use crate::{check_if_inbounds, check1d_data};

/// Linear Interpolator.
///
/// This interpolation method does not require any additional memory.
#[doc(alias = "gsl_interp_linear")]
#[derive(Debug, Clone)]
pub struct LinearInterpolator;

impl LinearInterpolator {
    /// The minimum required number of points.
    #[doc(alias = "gsl_interp_min_size")]
    const MIN_SIZE: usize = 2;

    /// Constructs a Linear Interpolator.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = LinearInterpolator::build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// - [`InterpolationError::UnsortedDataset`]: `xa` is not monotonically increasing.
    /// - [`InterpolationError::DatasetMismatch`]: `xa` and `ya` do not have the same length.
    /// - [`InterpolationError::NotEnoughPoints`]: length of `xa` is less that 2.
    #[doc(alias = "gsl_interp_init")]
    pub fn build(xa: &[f64], ya: &[f64]) -> Result<Self, InterpolationError> {
        check1d_data(xa, ya, Self::MIN_SIZE)?;
        Ok(Self)
    }
}

impl Interpolation for LinearInterpolator {
    fn eval(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        check_if_inbounds(xa, x)?;
        let index = acc.find(xa, x);

        let xlo = xa[index];
        let xhi = xa[index + 1];
        let ylo = ya[index];
        let yhi = ya[index + 1];

        let dx = xhi - xlo;

        debug_assert!(dx > 0.0);
        Ok(ylo + (x - xlo) / dx * (yhi - ylo))
    }

    fn eval_deriv(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        check_if_inbounds(xa, x)?;
        let index = acc.find(xa, x);

        let xlo = xa[index];
        let xhi = xa[index + 1];
        let ylo = ya[index];
        let yhi = ya[index + 1];

        let dx = xhi - xlo;
        let dy = yhi - ylo;

        debug_assert!(dx > 0.0);
        Ok(dy / dx)
    }

    /// Always returns `0`.
    fn eval_deriv2(
        &self,
        xa: &[f64],
        _: &[f64],
        x: f64,
        _: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        check_if_inbounds(xa, x)?;
        Ok(0.0)
    }

    fn eval_integ(
        &self,
        xa: &[f64],
        ya: &[f64],
        a: f64,
        b: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        check_if_inbounds(xa, a)?;
        check_if_inbounds(xa, b)?;
        let index_a = acc.find(xa, a);
        let index_b = acc.find(xa, b);

        let mut result = 0.0;

        for i in index_a..=index_b {
            let xlo = xa[i];
            let xhi = xa[i + 1];
            let ylo = ya[i];
            let yhi = ya[i + 1];

            let dx = xhi - xlo;
            let d = (yhi - ylo) / dx;

            if dx == 0.0 {
                continue;
            }

            if (i == index_a) | (i == index_b) {
                let x1 = if i == index_a { a } else { xlo };
                let x2 = if i == index_b { b } else { xhi };
                result += (x2 - x1) * (ylo + 0.5 * d * ((x2 - xlo) + (x1 - xlo)));
            } else {
                result += 0.5 * dx * (ylo + yhi);
            }
        }
        Ok(result)
    }
}
