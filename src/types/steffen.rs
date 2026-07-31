//! Definition Steffen Interpolator.

use crate::Accelerator;
use crate::{Domain1dError, InterpolatorError};
use crate::{Interpolation, Interpolator};
use crate::{check_if_inbounds, check1d_data};

const MIN_SIZE: usize = 3;

/// Steffen Interpolator
///
/// Steffen’s method guarantees the monotonicity of the interpolating function between the given
/// data points. Therefore, minima and maxima can only occur exactly at the data points, and there
/// can never be spurious oscillations between data points. The interpolated function is piecewise
/// cubic in each interval. The resulting curve and its first derivative are guaranteed to be
/// continuous, but the second derivative may be discontinuous.
#[doc(alias = "gsl_steffen_interp")]
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct SteffenInterpolator {
    a: Box<[f64]>,
    b: Box<[f64]>,
    c: Box<[f64]>,
    d: Box<[f64]>,
}

impl Interpolator for SteffenInterpolator {
    /// Constructs a Cubic Interpolator.
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = SteffenInterpolator::build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// - [`InterpolatorError::UnsortedDataset`]: `xa` is not monotonically increasing.
    /// - [`InterpolatorError::DatasetMismatch`]: `xa` and `ya` do not have the same length.
    /// - [`InterpolatorError::NotEnoughPoints`]: length of `xa` is less that 3.
    fn build(xa: &[f64], ya: &[f64]) -> Result<Self, InterpolatorError> {
        check1d_data(xa, ya, MIN_SIZE)?;
        let size = xa.len();

        // First assign the interval and slopes for the left boundary. We use the "simplest
        // possibility" method described in the paper in section 2.2.
        let h0 = xa[1] - xa[0];
        let s0 = (ya[1] - ya[0]) / h0;

        // Stores the calculated y' values.
        let mut y_prime = Vec::with_capacity(size);
        y_prime.push(s0);

        // Calculate all necessary s, h, p, y' variables form 1 to `size-2` (0 to `size-2`
        // inclusive)
        for i in 1..(size - 1) {
            // Equation 6 in paper
            let hi = xa[i + 1] - xa[i];
            let him1 = xa[i] - xa[i - 1];

            // Equation 7 in paper
            let si = (ya[i + 1] - ya[i]) / hi;
            let sim1 = (ya[i] - ya[i - 1]) / him1;

            // Equation 8 in paper
            let pi = (sim1 * hi + si * him1) / (him1 + hi);

            let min1 = si.abs().min(0.5 * pi.abs());
            let min2 = sim1.abs().min(min1);
            y_prime.push((steffen_copysign(1.0, sim1) + steffen_copysign(1.0, si)) * min2);
        }

        // We also need y' for the rightmost boundary; we use the 'simplest possibility' method
        // described in the paper in section 2.2
        //
        //y` = s_{n-1}`
        y_prime.push((ya[size - 1] - ya[size - 2]) / (xa[size - 1] - xa[size - 2]));

        // Now we can calculate all the coefficients for the whole range
        let mut a = Vec::with_capacity(size - 1);
        let mut b = Vec::with_capacity(size - 1);
        let mut c = Vec::with_capacity(size - 1);
        let mut d = Vec::with_capacity(size - 1);

        for i in 0..(size - 1) {
            let hi = xa[i + 1] - xa[i];
            let si = (ya[i + 1] - ya[i]) / hi;

            // These are from equations 2-5 in the paper
            a.push((y_prime[i] + y_prime[i + 1] - 2.0 * si) / hi.powi(2));
            b.push((3.0 * si - 2.0 * y_prime[i] - y_prime[i + 1]) / hi);
            c.push(y_prime[i]);
            d.push(ya[i]);
        }

        Ok(SteffenInterpolator {
            a: a.into_boxed_slice(),
            b: b.into_boxed_slice(),
            c: c.into_boxed_slice(),
            d: d.into_boxed_slice(),
        })
    }

    fn name(&self) -> &str {
        "Steffen"
    }

    fn min_size(&self) -> usize {
        MIN_SIZE
    }
}

// ===============================================================================================

impl Interpolation for SteffenInterpolator {
    fn eval(
        &self,
        xa: &[f64],
        _: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        check_if_inbounds(xa, x)?;
        let index = acc.find(xa, x);

        let xlo = xa[index];
        let delx = x - xlo;
        let a = self.a[index];
        let b = self.b[index];
        let c = self.c[index];
        let d = self.d[index];

        Ok(d + delx * (c + delx * (b + delx * a)))
    }

    fn eval_deriv(
        &self,
        xa: &[f64],
        _: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        check_if_inbounds(xa, x)?;
        let index = acc.find(xa, x);

        let xlo = xa[index];
        let delx = x - xlo;
        let a = self.a[index];
        let b = self.b[index];
        let c = self.c[index];

        Ok(c + delx * (2.0 * b + delx * 3.0 * a))
    }

    fn eval_deriv2(
        &self,
        xa: &[f64],
        _: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        check_if_inbounds(xa, x)?;
        let index = acc.find(xa, x);

        let xlo = xa[index];
        let delx = x - xlo;
        let a = self.a[index];
        let b = self.b[index];

        Ok(6.0 * delx * a + 2.0 * b)
    }

    fn eval_integ(
        &self,
        xa: &[f64],
        _: &[f64],
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

            let dx = xhi - xlo;

            // If two x points are the same
            if dx == 0.0 {
                continue;
            }

            let x1 = if i == index_a { a - xlo } else { 0.0 };
            let x2 = if i == index_b { b - xlo } else { xhi - xlo };

            let x12 = x1.powi(2);
            let x22 = x2.powi(2);

            result += 0.25 * self.a[i] * (x22.powi(2) - x12.powi(2));
            result += (1.0 / 3.0) * self.b[i] * (x22 * x2 - x12 * x1);
            result += 0.5 * self.c[i] * (x22 - x12);
            result += self.d[i] * (x2 - x1);
        }

        Ok(result)
    }
}

fn steffen_copysign(x: f64, y: f64) -> f64 {
    if (x.is_sign_negative() & y.is_sign_positive()) | (x.is_sign_positive() & y.is_sign_negative())
    {
        -x
    } else {
        x
    }
}
