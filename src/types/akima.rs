//! Definitions of `Akima` and `AkimaPeriodic` Interpolators.

use std::collections::VecDeque;

use crate::{Accelerator, BuildInterpolator, Domain1dError, Interpolation, InterpolationError};
use crate::{check_if_inbounds, check1d_data, integ_eval};

/// Akima Interpolator.
///
/// Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded
/// corner algorithm of Wodicka.
#[doc(alias = "gsl_akima_interp")]
#[derive(Debug, Clone)]
pub struct AkimaInterpolator {
    coef: Box<[(f64, f64, f64)]>,
}

impl BuildInterpolator for AkimaInterpolator {
    const MIN_SIZE: usize = 5;

    /// Constructs an Akima Interpolator.
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let interp = AkimaInterpolator::build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// - [`InterpolationError::UnsortedDataset`]: `xa` is not monotonically increasing.
    /// - [`InterpolationError::DatasetMismatch`]: `xa` and `ya` do not have the same length.
    /// - [`InterpolationError::NotEnoughPoints`]: length of `xa` is less that 5.
    fn build(xa: &[f64], ya: &[f64]) -> Result<Self, InterpolationError> {
        check1d_data(xa, ya, Self::MIN_SIZE)?;

        let size = xa.len();

        // All m indices are shifted by +2
        let mut m = VecDeque::with_capacity(size);
        for i in 0..=size - 2 {
            m.push_back((ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]));
        }

        // Non-periodic boundary conditions
        m.push_front(2.0 * m[0] - m[1]);
        m.push_front(3.0 * m[1] - 2.0 * m[2]);
        m.push_back(2.0 * m[size] - m[size - 1]);
        m.push_back(3.0 * m[size] - 2.0 * m[size - 1]);
        let m = m.make_contiguous();

        let (b, c, d) = akima_calc(xa, m);

        let mut coef: Vec<(f64, f64, f64)> = Vec::with_capacity(size - 1);
        for i in 0..size - 1 {
            coef.push((b[i], c[i], d[i]));
        }

        Ok(AkimaInterpolator {
            coef: coef.into_boxed_slice(),
        })
    }
}

// ===============================================================================================

impl Interpolation for AkimaInterpolator {
    fn eval(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        akima_eval(xa, ya, &self.coef, x, acc)
    }

    fn eval_deriv(
        &self,
        xa: &[f64],
        _: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        akima_eval_deriv(xa, &self.coef, x, acc)
    }

    fn eval_deriv2(
        &self,
        xa: &[f64],
        _: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        akima_eval_deriv2(xa, &self.coef, x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[f64],
        ya: &[f64],
        a: f64,
        b: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        akima_eval_integ(xa, ya, &self.coef, a, b, acc)
    }
}

//=================================================================================================

/// Akima Interpolator.
///
/// Non-rounded Akima spline with natural boundary conditions. This method uses the non-rounded
/// corner algorithm of Wodicka.
#[doc(alias = "gsl_interp_akima_periodic")]
#[derive(Debug, Clone)]
pub struct AkimaPeriodicInterpolator {
    coef: Box<[(f64, f64, f64)]>,
}

impl BuildInterpolator for AkimaPeriodicInterpolator {
    const MIN_SIZE: usize = 5;

    /// Constructs an Akima Periodic Interpolator.
    ///
    /// ## Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    /// let interp = AkimaPeriodicInterpolator::build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// - [`InterpolationError::UnsortedDataset`]: `xa` is not monotonically increasing.
    /// - [`InterpolationError::DatasetMismatch`]: `xa` and `ya` do not have the same length.
    /// - [`InterpolationError::NotEnoughPoints`]: length of `xa` is less that 5.
    fn build(xa: &[f64], ya: &[f64]) -> Result<Self, InterpolationError> {
        check1d_data(xa, ya, Self::MIN_SIZE)?;

        let size = xa.len();

        // All m indices are shifted by +2
        let mut m = VecDeque::with_capacity(size + 3);
        for i in 0..=size - 2 {
            m.push_back((ya[i + 1] - ya[i]) / (xa[i + 1] - xa[i]));
        }

        // Periodic boundary conditions
        m.push_front(m[size - 1 - 1]);
        m.push_front(m[size - 1 - 1]);
        m.push_back(m[2]);
        m.push_back(m[3]);
        let m = m.make_contiguous();

        let (b, c, d) = akima_calc(xa, m);

        let mut coef: Vec<(f64, f64, f64)> = Vec::with_capacity(size - 1);
        for i in 0..size - 1 {
            coef.push((b[i], c[i], d[i]));
        }

        Ok(AkimaPeriodicInterpolator {
            coef: coef.into_boxed_slice(),
        })
    }
}

// ===============================================================================================

impl Interpolation for AkimaPeriodicInterpolator {
    fn eval(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        akima_eval(xa, ya, &self.coef, x, acc)
    }

    fn eval_deriv(
        &self,
        xa: &[f64],
        _: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        akima_eval_deriv(xa, &self.coef, x, acc)
    }

    fn eval_deriv2(
        &self,
        xa: &[f64],
        _: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        akima_eval_deriv2(xa, &self.coef, x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[f64],
        ya: &[f64],
        a: f64,
        b: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        akima_eval_integ(xa, ya, &self.coef, a, b, acc)
    }
}

//=================================================================================================

fn akima_eval(
    xa: &[f64],
    ya: &[f64],
    coef: &[(f64, f64, f64)],
    x: f64,
    acc: &mut Accelerator,
) -> Result<f64, Domain1dError> {
    check_if_inbounds(xa, x)?;
    let index = acc.find(xa, x);

    let xlo = xa[index];
    let delx = x - xlo;
    let (b, c, d) = coef[index];

    Ok(ya[index] + delx * (b + delx * (c + d * delx)))
}

fn akima_eval_deriv(
    xa: &[f64],
    coef: &[(f64, f64, f64)],
    x: f64,
    acc: &mut Accelerator,
) -> Result<f64, Domain1dError> {
    check_if_inbounds(xa, x)?;

    let index = acc.find(xa, x);

    let xlo = xa[index];
    let delx = x - xlo;
    let (b, c, d) = coef[index];

    Ok(b + delx * (2.0 * c + 3.0 * d * delx))
}

fn akima_eval_deriv2(
    xa: &[f64],
    coef: &[(f64, f64, f64)],
    x: f64,
    acc: &mut Accelerator,
) -> Result<f64, Domain1dError> {
    check_if_inbounds(xa, x)?;

    let index = acc.find(xa, x);

    let xlo = xa[index];
    let delx = x - xlo;
    let (_, c, d) = coef[index];

    Ok(2.0 * c + 6.0 * d * delx)
}

fn akima_eval_integ(
    xa: &[f64],
    ya: &[f64],
    coef: &[(f64, f64, f64)],
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

        let dx = xhi - xlo;

        // If two x points are the same
        if dx == 0.0 {
            continue;
        }

        let (bs, cs, ds) = coef[i];
        if (i == index_a) | (i == index_b) {
            let x1 = if i == index_a { a } else { xlo };
            let x2 = if i == index_b { b } else { xhi };
            result += integ_eval(ylo, bs, cs, ds, xlo, x1, x2);
        } else {
            result += dx * (ylo + dx * (0.5 * bs + dx * ((1.0 / 3.0) * cs + 0.25 * ds * dx)))
        }
    }
    Ok(result)
}

/// Common Calculation.
#[expect(clippy::type_complexity, reason = "only used here")]
fn akima_calc(xa: &[f64], m: &[f64]) -> (Box<[f64]>, Box<[f64]>, Box<[f64]>) {
    let size = xa.len();
    let mut b = Vec::with_capacity(size - 1);
    let mut c = Vec::with_capacity(size - 1);
    let mut d = Vec::with_capacity(size - 1);

    for i in 0..size - 1 {
        let ne = (m[i + 3] - m[i + 2]).abs() + (m[i + 1] - m[i]).abs();
        if ne == 0.0 {
            b.push(m[i + 2]);
            c.push(0.0);
            d.push(0.0);
        } else {
            let hi = xa[i + 1] - xa[i];
            let nenext = (m[i + 4] - m[i + 3]).abs() + (m[i + 2] - m[i + 1]).abs();
            let ai = (m[i + 1] - m[i]).abs() / ne;
            let ai_plus1: f64;
            let tli_plus1 = if nenext == 0.0 {
                m[i + 2]
            } else {
                ai_plus1 = (m[i + 2] - m[i + 1]).abs() / nenext;
                (1.0 - ai_plus1) * m[i + 2] + ai_plus1 * m[i + 3]
            };
            b.push((1.0 - ai) * m[i + 1] + ai * m[i + 2]);
            c.push((3.0 * m[i + 2] - 2.0 * b[i] - tli_plus1) / hi);
            d.push((b[i] + tli_plus1 - 2.0 * m[i + 2]) / hi.powi(2));
        }
    }
    (
        b.into_boxed_slice(),
        c.into_boxed_slice(),
        d.into_boxed_slice(),
    )
}
