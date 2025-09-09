use crate::Accelerator;
use crate::DomainError;
use crate::InterpType;
use crate::Interpolation;
use crate::InterpolationError;
use crate::types::utils::{check_if_inbounds, check1d_data};

const MIN_SIZE: usize = 3;

/// Steffen Interpolation type.
///
/// Steffen’s method guarantees the monotonicity of the interpolating function between the given
/// data points. Therefore, minima and maxima can only occur exactly at the data points, and there
/// can never be spurious oscillations between data points. The interpolated function is piecewise
/// cubic in each interval. The resulting curve and its first derivative are guaranteed to be
/// continuous, but the second derivative may be discontinuous.
pub struct Steffen;

impl<T> InterpType<T> for Steffen
where
    T: crate::Num,
{
    type Interpolation = SteffenInterp<T>;

    /// Constructs a Cubic Interpolator.
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = Steffen.build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    fn build(&self, xa: &[T], ya: &[T]) -> Result<SteffenInterp<T>, InterpolationError> {
        check1d_data(xa, ya, MIN_SIZE)?;
        let size = xa.len();

        // First assign the interval and slopes for the left boundary. We use the "simplest
        // possibility" method described in the paper in section 2.2.
        let h0 = xa[1] - xa[0];
        let s0 = (ya[1] - ya[0]) / h0;

        // Stores the calculated y' values.
        let mut y_prime = Vec::<T>::with_capacity(size);
        y_prime.push(s0);

        // Calculate all necessary s, h, p, y' variables form 1 to `size-2` (0 to `size-2`
        // inclusive)
        let one = T::one();
        let half = T::from(0.5).unwrap();

        for i in 1..(size - 1) {
            // Equation 6 in paper
            let hi = xa[i + 1] - xa[i];
            let him1 = xa[i] - xa[i - 1];

            // Equation 7 in paper
            let si = (ya[i + 1] - ya[i]) / hi;
            let sim1 = (ya[i] - ya[i - 1]) / him1;

            // Equation 8 in paper
            let pi = (sim1 * hi + si * him1) / (him1 + hi);

            let min1 = si.abs().min(half * pi.abs());
            let min2 = sim1.abs().min(min1);
            y_prime.push((steffen_copysign(one, sim1) + steffen_copysign(one, si)) * min2);
        }

        // We also need y' for the rightmost boundary; we use the 'simplest possibility' method
        // described in the paper in section 2.2
        //
        //y` = s_{n-1}`
        y_prime.push((ya[size - 1] - ya[size - 2]) / (xa[size - 1] - xa[size - 2]));

        // Now we can calculate all the coefficients for the whole range
        let mut a = Vec::<T>::with_capacity(size - 1);
        let mut b = Vec::<T>::with_capacity(size - 1);
        let mut c = Vec::<T>::with_capacity(size - 1);
        let mut d = Vec::<T>::with_capacity(size - 1);

        let two = T::from(2.0).unwrap();
        let three = T::from(3.0).unwrap();
        for i in 0..(size - 1) {
            let hi = xa[i + 1] - xa[i];
            let si = (ya[i + 1] - ya[i]) / hi;

            // These are from equations 2-5 in the paper
            a.push((y_prime[i] + y_prime[i + 1] - two * si) / hi.powi(2));
            b.push((three * si - two * y_prime[i] - y_prime[i + 1]) / hi);
            c.push(y_prime[i]);
            d.push(ya[i]);
        }

        let state = SteffenInterp {
            a,
            b,
            c,
            d,
            y_prime,
        };

        Ok(state)
    }

    fn name(&self) -> &str {
        "Steffen"
    }

    fn min_size(&self) -> usize {
        MIN_SIZE
    }
}

// ===============================================================================================

/// Steffen Interpolator.
///
/// Provides all the evaluation methods.
///
/// Should be constructed through the [`Steffen`] type.
#[allow(dead_code)]
pub struct SteffenInterp<T>
where
    T: crate::Num,
{
    a: Vec<T>,
    b: Vec<T>,
    c: Vec<T>,
    d: Vec<T>,
    y_prime: Vec<T>,
}

impl<T> Interpolation<T> for SteffenInterp<T>
where
    T: crate::Num,
{
    #[allow(unused_variables)]
    fn eval(&self, xa: &[T], ya: &[T], x: T, acc: &mut Accelerator) -> Result<T, DomainError> {
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

    #[allow(unused_variables)]
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
        let delx = x - xlo;
        let a = self.a[index];
        let b = self.b[index];
        let c = self.c[index];

        let two = T::from(2.0).unwrap();
        let three = T::from(3.0).unwrap();

        Ok(c + delx * (two * b + delx * three * a))
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
        let index = acc.find(xa, x);

        let xlo = xa[index];
        let delx = x - xlo;
        let a = self.a[index];
        let b = self.b[index];

        let two = T::from(2.0).unwrap();
        let six = T::from(6.0).unwrap();

        Ok(six * delx * a + two * b)
    }

    #[allow(unused_variables)]
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

        let quarter = T::from(0.25).unwrap();
        let half = T::from(0.5).unwrap();
        let third = T::from(1.0 / 3.0).unwrap();
        let mut result = T::zero();

        for i in index_a..=index_b {
            let xlo = xa[i];
            let xhi = xa[i + 1];

            let dx = xhi - xlo;

            // If two x points are the same
            if dx.is_zero() {
                continue;
            }

            let x1 = if i == index_a { a - xlo } else { T::zero() };
            let x2 = if i == index_b { b - xlo } else { xhi - xlo };

            let x12 = x1.powi(2);
            let x22 = x2.powi(2);

            result += quarter * self.a[i] * (x22.powi(2) - x12.powi(2));
            result += third * self.b[i] * (x22 * x2 - x12 * x1);
            result += half * self.c[i] * (x22 - x12);
            result += self.d[i] * (x2 - x1);
        }

        Ok(result)
    }
}

fn steffen_copysign<T>(x: T, y: T) -> T
where
    T: crate::Num,
{
    if (x.is_sign_negative() & y.is_sign_positive()) | (x.is_sign_positive() & y.is_sign_negative())
    {
        -x
    } else {
        x
    }
}
