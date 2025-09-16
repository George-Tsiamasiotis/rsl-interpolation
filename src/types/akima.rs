use std::collections::VecDeque;

use crate::Accelerator;
use crate::DomainError;
use crate::InterpType;
use crate::Interpolation;
use crate::InterpolationError;
use crate::types::utils::integ_eval;
use crate::types::utils::{check_if_inbounds, check1d_data};

const MIN_SIZE: usize = 5;

/// Akima Interpolation type.
///
/// Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded
/// corner algorithm of Wodicka.
#[doc(alias = "gsl_interp_akima")]
pub struct Akima;

impl<T> InterpType<T> for Akima
where
    T: crate::Num,
{
    type Interpolation = AkimaInterp<T>;

    /// Constructs an Akima Interpolator.
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let interp = Akima.build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    fn build(&self, xa: &[T], ya: &[T]) -> Result<AkimaInterp<T>, InterpolationError> {
        check1d_data(xa, ya, MIN_SIZE)?;

        let size = xa.len();
        let two = T::from(2.0).unwrap();
        let three = T::from(3.0).unwrap();

        // All m indices are shifted by +2
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

        let state = AkimaInterp { b, c, d, m };
        Ok(state)
    }

    fn name(&self) -> &str {
        "Akima"
    }

    fn min_size(&self) -> usize {
        MIN_SIZE
    }
}

// ===============================================================================================

/// Akima Interpolator.
///
/// Provides all the evaluation methods.
///
/// Should be constructed through the [`Akima`] type.
#[allow(dead_code)]
#[doc(alias = "gsl_akima_interp")]
pub struct AkimaInterp<T>
where
    T: crate::Num,
{
    b: Vec<T>,
    c: Vec<T>,
    d: Vec<T>,
    m: Vec<T>,
}

impl<T> Interpolation<T> for AkimaInterp<T>
where
    T: crate::Num,
{
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

/// Akima Periodic Interpolation type.
///
/// Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded
/// corner algorithm of Wodicka.
pub struct AkimaPeriodic;

impl<T> InterpType<T> for AkimaPeriodic
where
    T: crate::Num,
{
    type Interpolation = AkimaPeriodicInterp<T>;

    /// Constructs an Akima Periodic Interpolator.
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let interp = AkimaPeriodic.build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    fn build(&self, xa: &[T], ya: &[T]) -> Result<AkimaPeriodicInterp<T>, InterpolationError> {
        check1d_data(xa, ya, MIN_SIZE)?;

        let size = xa.len();

        // All m indices are shifted by +2
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

        let state = AkimaPeriodicInterp { b, c, d, m };
        Ok(state)
    }

    fn name(&self) -> &str {
        "Akima Periodic"
    }

    fn min_size(&self) -> usize {
        MIN_SIZE
    }
}

// ===============================================================================================

/// Akima Interpolator.
///
/// Provides all the evaluation methods.
///
/// Should be constructed through the [`Akima`] type.
#[allow(dead_code)]
#[doc(alias = "gsl_interp_akima_periodic")]
pub struct AkimaPeriodicInterp<T>
where
    T: crate::Num,
{
    c: Vec<T>,
    b: Vec<T>,
    d: Vec<T>,
    m: Vec<T>,
}

impl<T> Interpolation<T> for AkimaPeriodicInterp<T>
where
    T: crate::Num,
{
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
    T: crate::Num,
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
    T: crate::Num,
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

#[inline(always)]
fn akima_eval_deriv2<T>(
    xa: &[T],
    state: (&[T], &[T]),
    x: T,
    acc: &mut Accelerator,
) -> Result<T, DomainError>
where
    T: crate::Num,
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
    T: crate::Num,
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
            result += integ_eval(ylo, bs[i], cs[i], ds[i], xlo, x1, x2);
        } else {
            result += dx * (ylo + dx * (half * bs[i] + dx * (third * cs[i] + quarter * ds[i] * dx)))
        }
    }
    Ok(result)
}

/// Common Calculation
fn akima_calc<T>(xa: &[T], m: &[T]) -> (Vec<T>, Vec<T>, Vec<T>)
where
    T: crate::Num,
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
