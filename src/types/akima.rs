//! Implementor for Akima Interpolator.

use std::collections::VecDeque;
use std::fmt::Debug;

use crate::Accelerator;
use crate::DomainError;
use crate::Interpolation;
use crate::types::utils::integ_eval;
use crate::types::utils::{check_if_inbounds, check1d_data};

/// Akima Spline.
///
/// Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded
/// corner algorithm of Wodicka.
///
/// ## Example
///
/// ```
/// # use rsl_interpolation::Interpolation;
/// # use rsl_interpolation::InterpolationError;
/// # use rsl_interpolation::Akima;
/// # use rsl_interpolation::Accelerator;
/// #
/// # fn main() -> Result<(), InterpolationError>{
/// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
/// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
/// let interp = Akima::new(&xa, &ya)?;
/// # Ok(())
/// # }
/// ```
#[allow(dead_code)]
pub struct Akima<T>
where
    T: num::Float + std::fmt::Debug,
{
    b: Vec<T>,
    c: Vec<T>,
    d: Vec<T>,
    m: Vec<T>,
}

impl<T> Interpolation<T> for Akima<T>
where
    T: num::Float + Debug,
{
    const MIN_SIZE: usize = 5;
    const NAME: &'static str = "akima";

    fn new(xa: &[T], ya: &[T]) -> Result<Self, crate::InterpolationError>
    where
        Self: Sized,
    {
        check1d_data(xa, ya, Self::MIN_SIZE)?;

        let size = xa.len();
        let two = T::from(2.0).unwrap();
        let three = T::from(3.0).unwrap();

        // All m indeces are shifted by +2
        let mut m = VecDeque::<T>::with_capacity(size);
        for i in 0..=size - 2 {
            m.push_back((ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]));
        }

        // Non-periodic boundary conditions
        m.push_front(two * m[0] - m[1]);
        m.push_front(three * m[1] - two * m[2]);
        m.push_back(two * m[size] - m[size - 1]);
        m.push_back(three * m[size] - two * m[size - 1]);
        let m = m.make_contiguous().to_vec();

        let (b, c, d) = akima_calc(xa, &m);

        let akima = Akima { b, c, d, m };
        Ok(akima)
    }

    fn eval(&self, xa: &[T], ya: &[T], x: T, acc: &mut Accelerator) -> Result<T, DomainError> {
        akima_eval(xa, ya, (&self.b, &self.c, &self.d), x, acc)
    }

    #[allow(unused_variables)]
    fn eval_deriv(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        akima_eval_deriv(xa, (&self.b, &self.c, &self.d), x, acc)
    }

    #[allow(unused_variables)]
    fn eval_deriv2(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        akima_eval_deriv2(xa, (&self.c, &self.d), x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[T],
        ya: &[T],
        a: T,
        b: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        akima_eval_integ(xa, ya, (&self.b, &self.c, &self.d), a, b, acc)
    }
}

//=================================================================================================

/// Akima Periodic Spline.
///
/// Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded
/// corner algorithm of Wodicka.
///
/// ## Example
///
/// ```
/// # use rsl_interpolation::Interpolation;
/// # use rsl_interpolation::InterpolationError;
/// # use rsl_interpolation::AkimaPeriodic;
/// # use rsl_interpolation::Accelerator;
/// #
/// # fn main() -> Result<(), InterpolationError>{
/// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
/// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
/// let interp = AkimaPeriodic::new(&xa, &ya)?;
/// # Ok(())
/// # }
/// ```
#[allow(dead_code)]
pub struct AkimaPeriodic<T>
where
    T: num::Float + Debug,
{
    c: Vec<T>,
    b: Vec<T>,
    d: Vec<T>,
    m: Vec<T>,
}

impl<T> Interpolation<T> for AkimaPeriodic<T>
where
    T: num::Float + Debug,
{
    const MIN_SIZE: usize = 5;
    const NAME: &'static str = "akima-periodic";

    fn new(xa: &[T], ya: &[T]) -> Result<Self, crate::InterpolationError>
    where
        Self: Sized,
    {
        check1d_data(xa, ya, Self::MIN_SIZE)?;

        let size = xa.len();

        // All m indeces are shifted by +2
        let mut m = VecDeque::<T>::with_capacity(size + 3);
        for i in 0..=size - 2 {
            m.push_back((ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]));
        }

        // Periodic boundary conditions
        m.push_front(m[size - 1 - 1]);
        m.push_front(m[size - 1 - 1]);
        m.push_back(m[2]);
        m.push_back(m[3]);
        let m = m.make_contiguous().to_vec();

        let (b, c, d) = akima_calc(xa, &m);

        let akima_periodic = AkimaPeriodic { b, c, d, m };
        Ok(akima_periodic)
    }

    fn eval(&self, xa: &[T], ya: &[T], x: T, acc: &mut Accelerator) -> Result<T, DomainError> {
        akima_eval(xa, ya, (&self.b, &self.c, &self.d), x, acc)
    }

    #[allow(unused_variables)]
    fn eval_deriv(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        akima_eval_deriv(xa, (&self.b, &self.c, &self.d), x, acc)
    }

    #[allow(unused_variables)]
    fn eval_deriv2(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        akima_eval_deriv2(xa, (&self.c, &self.d), x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[T],
        ya: &[T],
        a: T,
        b: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        akima_eval_integ(xa, ya, (&self.b, &self.c, &self.d), a, b, acc)
    }
}

//=================================================================================================

fn akima_eval<T>(
    xa: &[T],
    ya: &[T],
    state: (&[T], &[T], &[T]),
    x: T,
    acc: &mut Accelerator,
) -> Result<T, DomainError>
where
    T: num::Float + Debug,
{
    check_if_inbounds(xa, x)?;
    let index = acc.find(xa, x);

    let xlo = xa[index];
    let delx = x - xlo;
    let b = state.0[index];
    let c = state.1[index];
    let d = state.2[index];

    Ok(ya[index] + delx * (b + delx * (c + d * delx)))
}

fn akima_eval_deriv<T>(
    xa: &[T],
    state: (&[T], &[T], &[T]),
    x: T,
    acc: &mut Accelerator,
) -> Result<T, DomainError>
where
    T: num::Float + Debug,
{
    check_if_inbounds(xa, x)?;
    let two = T::from(2).unwrap();
    let three = T::from(3).unwrap();

    let index = acc.find(xa, x);

    let xlo = xa[index];
    let delx = x - xlo;
    let b = state.0[index];
    let c = state.1[index];
    let d = state.2[index];

    Ok(b + delx * (two * c + three * d * delx))
}

fn akima_eval_deriv2<T>(
    xa: &[T],
    state: (&[T], &[T]),
    x: T,
    acc: &mut Accelerator,
) -> Result<T, DomainError>
where
    T: num::Float + Debug,
{
    check_if_inbounds(xa, x)?;
    let two = T::from(2).unwrap();
    let six = T::from(6).unwrap();

    let index = acc.find(xa, x);

    let xlo = xa[index];
    let delx = x - xlo;
    let c = state.0[index];
    let d = state.1[index];

    Ok(two * c + six * d * delx)
}

fn akima_eval_integ<T>(
    xa: &[T],
    ya: &[T],
    state: (&[T], &[T], &[T]),
    a: T,
    b: T,
    acc: &mut Accelerator,
) -> Result<T, DomainError>
where
    T: num::Float + Debug,
{
    check_if_inbounds(xa, a)?;
    check_if_inbounds(xa, b)?;
    let index_a = acc.find(xa, a);
    let index_b = acc.find(xa, b);
    let bs = state.0;
    let cs = state.1;
    let ds = state.2;

    let quarter = T::from(0.25).unwrap();
    let half = T::from(0.5).unwrap();
    let third = T::from(1.0 / 3.0).unwrap();
    let mut result = T::zero();

    for i in index_a..=index_b {
        let xlo = xa[i];
        let xhi = xa[i + 1];
        let ylo = ya[i];

        let dx = xhi - xlo;

        // If two x points are the same
        if dx.is_zero() {
            continue;
        }

        if (i == index_a) | (i == index_b) {
            let x1 = if i == index_a { a } else { xlo };
            let x2 = if i == index_b { b } else { xhi };
            result = result + integ_eval(ylo, bs[i], cs[i], ds[i], xlo, x1, x2);
        } else {
            result = result
                + dx * (ylo + dx * (half * bs[i] + dx * (third * cs[i] + quarter * ds[i] * dx)))
        }
    }
    Ok(result)
}

/// Common Calculation
fn akima_calc<T>(xa: &[T], m: &[T]) -> (Vec<T>, Vec<T>, Vec<T>)
where
    T: num::Float + Debug,
{
    let size = xa.len();
    let two = T::from(2.0).unwrap();
    let three = T::from(3.0).unwrap();
    let mut b = Vec::<T>::with_capacity(size - 1);
    let mut c = Vec::<T>::with_capacity(size - 1);
    let mut d = Vec::<T>::with_capacity(size - 1);

    for i in 0..size - 1 {
        let ne = (m[i + 3] - m[i + 2]).abs() + (m[i + 1] - m[i]).abs();
        if ne.is_zero() {
            b.push(m[i + 2]);
            c.push(T::zero());
            d.push(T::zero());
        } else {
            let hi = xa[i + 1] - xa[i];
            let nenext = (m[i + 4] - m[i + 3]).abs() + (m[i + 2] - m[i + 1]).abs();
            let ai = (m[i + 1] - m[i]).abs() / ne;
            let ai_plus1: T;
            let tli_plus1: T;
            if nenext.is_zero() {
                tli_plus1 = m[i + 2];
            } else {
                ai_plus1 = (m[i + 2] - m[i + 1]).abs() / nenext;
                tli_plus1 = (T::one() - ai_plus1) * m[i + 2] + ai_plus1 * m[i + 3];
            }
            b.push((T::one() - ai) * m[i + 1] + ai * m[i + 2]);
            c.push((three * m[i + 2] - two * b[i] - tli_plus1) / hi);
            d.push((b[i] + tli_plus1 - two * m[i + 2]) / hi.powi(2));
        }
    }
    (b, c, d)
}
