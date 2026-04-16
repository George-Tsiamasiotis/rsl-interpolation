//! Definitions of `Cubic` and `CubicPeriodic` Interpolators.

use ndarray::Array1;
use ndarray_linalg::{MatrixLayout, SolveTridiagonal, Tridiagonal};

use crate::{Accelerator, BuildInterpolator, Domain1dError, Interpolation, InterpolationError};
use crate::{check_if_inbounds, check1d_data, diff, integ_eval};

/// Cubic Interpolator.
///
/// Cubic Interpolation with natural boundary conditions. The resulting curve is piecewise cubic on each
/// interval, with matching first and second derivatives at the supplied data-points. The second
/// derivative is chosen to be zero at the first and last point.
///
/// ## Reference
///
/// Numerical Algorithms with C - Gisela Engeln-Mullges, Frank Uhlig - 1996 -
/// Algorithm 10.1, pg 254.
#[doc(alias = "gsl_interp_cspline")]
#[derive(Debug, Clone)]
pub struct CubicInterpolator {
    c: Box<[f64]>,
}

impl BuildInterpolator for CubicInterpolator {
    const MIN_SIZE: usize = 3;

    /// Constructs a Cubic Interpolator.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicInterpolator::build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// - [`InterpolationError::UnsortedDataset`]: `xa` is not monotonically increasing.
    /// - [`InterpolationError::DatasetMismatch`]: `xa` and `ya` do not have the same length.
    /// - [`InterpolationError::NotEnoughPoints`]: length of `xa` is less that 3.
    /// - [`InterpolationError::BLASTridiagError`]: Error when solving the tridiagonal system.
    fn build(xa: &[f64], ya: &[f64]) -> Result<Self, InterpolationError> {
        check1d_data(xa, ya, Self::MIN_SIZE)?;

        // Engeln-Mullges G. - Uhlig F.: Algorithm 10.1, pg 254
        let sys_size = xa.len() - 2;

        let h = diff(xa);
        debug_assert_eq!(h.len(), xa.len() - 1);

        // Ac=g setup
        let mut g = Vec::with_capacity(sys_size);
        let mut diag = Vec::with_capacity(sys_size);
        let mut offdiag = Vec::with_capacity(sys_size);
        for i in 0..sys_size {
            g.push(if h[i] == 0.0 {
                0.0
            } else {
                3.0 * (ya[i + 2] - ya[i + 1]) / h[i + 1] - 3.0 * (ya[i + 1] - ya[i]) / h[i]
            });
            diag.push(2.0 * (h[i] + h[i + 1]));
            offdiag.push(h[i + 1]);
        }
        // The last element of offdiag is not actually valid, by definition. Popping it is not
        // really needed though, since the solver ignores it. However, it is needed in the
        // CubicPeriodic case, since it represents the cyclical term.
        offdiag.pop();
        debug_assert_eq!(diag.len(), offdiag.len() + 1);

        let matrix = Tridiagonal {
            l: MatrixLayout::C {
                row: (sys_size) as i32,
                lda: (sys_size) as i32,
            },
            d: diag.clone(),
            dl: offdiag.clone(),
            du: offdiag.clone(),
        };

        // Ac=g solving
        let mut c = Vec::with_capacity(xa.len());
        c.push(0.0);
        if sys_size == 1 {
            c.push(g[0] / diag[0]);
        } else {
            let coeffs = match matrix.solve_tridiagonal(&Array1::from_vec(g.clone())) {
                Ok(coeffs) => coeffs,
                Err(err) => {
                    return Err(InterpolationError::BLASTridiagError {
                        which_interp: "Cubic".into(),
                        source: err,
                    });
                }
            };
            c = [c, coeffs.to_vec()].concat();
        }
        c.push(0.0);

        // g, diag, and offdiag are only needed for the calculation of c and are not used again
        Ok(CubicInterpolator {
            c: c.into_boxed_slice(),
        })
    }
}

// ===============================================================================================

impl Interpolation for CubicInterpolator {
    fn eval(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        cubic_eval(xa, ya, &self.c, x, acc)
    }

    fn eval_deriv(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        cubic_eval_deriv(xa, ya, &self.c, x, acc)
    }

    fn eval_deriv2(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        cubic_eval_deriv2(xa, ya, &self.c, x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[f64],
        ya: &[f64],
        a: f64,
        b: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        cubic_eval_integ(xa, ya, &self.c, a, b, acc)
    }
}

//=================================================================================================

/// Cubic Periodic interpolator.
///
/// Cubic Spline with periodic boundary conditions. The resulting curve is piecewise cubic on each
/// interval, with matching first and second derivatives at the supplied data-points. The
/// derivatives at the first and last points are also matched. Note that the last point in the data
/// must have the same y-value as the first point, otherwise the resulting periodic interpolation
/// will have a discontinuity at the boundary.
#[doc(alias = "gsl_interp_cspline_periodic")]
#[derive(Debug, Clone)]
pub struct CubicPeriodicInterpolator {
    c: Box<[f64]>,
}

impl BuildInterpolator for CubicPeriodicInterpolator {
    const MIN_SIZE: usize = 3;

    /// Constructs a Cubic Periodic Interpolator.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// let interp = CubicPeriodicInterpolator::build(&xa, &ya)?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Panics
    ///
    /// This type currently panics if built with more than 3 points, since it requires a
    /// cyclically tridiagonal matrix solver which is currently not implemented by
    /// `ndarray_linalg`.
    ///
    /// # Errors
    ///
    /// - [`InterpolationError::UnsortedDataset`]: `xa` is not monotonically increasing.
    /// - [`InterpolationError::DatasetMismatch`]: `xa` and `ya` do not have the same length.
    /// - [`InterpolationError::NotEnoughPoints`]: length of `xa` is less that 3.
    /// - [`InterpolationError::BLASTridiagError`]: Error when solving the tridiagonal system.
    #[expect(clippy::panic_in_result_fn, reason = "unimplemented solver")]
    fn build(xa: &[f64], ya: &[f64]) -> Result<Self, InterpolationError> {
        check1d_data(xa, ya, Self::MIN_SIZE)?;

        // Engeln-Mullges G. - Uhlig F.: Algorithm 10.2, pg 255
        let sys_size = xa.len() - 1;

        let h = diff(xa);
        debug_assert!(h.len() == xa.len() - 1);

        // Ac=g setup
        let mut c = Vec::with_capacity(xa.len());
        let mut g = Vec::with_capacity(sys_size);
        let mut diag = Vec::with_capacity(sys_size);
        let mut offdiag = Vec::with_capacity(sys_size);

        if sys_size == 2 {
            let h0 = xa[1] - xa[0];
            let h1 = xa[2] - xa[1];

            let a = 2.0 * (h0 + h1);
            let b = h0 + h1;

            g.push(3.0 * ((ya[2] - ya[1]) / h1 - (ya[1] - ya[0]) / h0));
            g.push(3.0 * ((ya[1] - ya[2]) / h0 - (ya[2] - ya[1]) / h1));

            let det = 3.0 * (h0 + h1) * (h0 + h1);
            c.push((-b * g[0] + a * g[1]) / det);
            c.push((a * g[0] - b * g[1]) / det);
            c.push(c[0]);
        } else {
            // Same as in Cubic case
            for i in 0..sys_size - 1 {
                g.push(if h[i] == 0.0 {
                    0.0
                } else {
                    3.0 * (ya[i + 2] - ya[i + 1]) / h[i + 1] - 3.0 * (ya[i + 1] - ya[i]) / h[i]
                });
                diag.push(2.0 * (h[i] + h[i + 1]));
                offdiag.push(h[i + 1]);
            }

            // But we must add the last point
            let i = sys_size - 1;
            let hi = xa[i + 1] - xa[i];
            let hiplus1 = xa[1] - xa[0];
            let ydiffi = ya[i + 1] - ya[i];
            let ydiffplus1 = ya[1] - ya[0];
            let gi = if !(hi == 0.0) { 1.0 / hi } else { 0.0 };
            let giplus1 = if !(hiplus1 == 0.0) {
                1.0 / hiplus1
            } else {
                0.0
            };
            offdiag.push(hiplus1);
            diag.push(2.0 * (hiplus1 + hi));
            g.push(3.0 * (ydiffplus1 * giplus1 - ydiffi * gi));
            // offdiag's last element represents the cyclical term
            debug_assert_eq!(diag.len(), offdiag.len());

            let matrix = Tridiagonal {
                l: MatrixLayout::C {
                    row: (sys_size) as i32,
                    lda: (sys_size) as i32,
                },
                d: diag.clone(),
                dl: offdiag.clone(),
                du: offdiag.clone(),
            };

            // Ac=g solving
            c.push(0.0);
            if sys_size == 1 {
                c.push(g[0] / diag[0]);
            } else {
                // This must solve a cyclically tridiagonal matrix, but its not implemented yet :(
                // The corner element is stored at the end of the offdiag vec.
                let coeffs = match matrix.solve_tridiagonal(&Array1::from_vec(g.clone())) {
                    Ok(coeffs) => coeffs,
                    Err(err) => {
                        return Err(InterpolationError::BLASTridiagError {
                            which_interp: "Cubic Periodic".into(),
                            source: err,
                        });
                    }
                };
                c = [c, coeffs.to_vec()].concat();
            }
            c[0] = c[sys_size];
            panic!(
                "\nNot implemented: Cubic Periodic Splines with more than 3 points require a solver for\
                cyclically tridiagonal matrices, which is currently not implemented by ndarray_linalg.\n"
            )
        }

        // g, diag, and offdiag are only needed for the calculation of c and are not used again
        Ok(CubicPeriodicInterpolator {
            c: c.into_boxed_slice(),
        })
    }
}

// ===============================================================================================

impl Interpolation for CubicPeriodicInterpolator {
    fn eval(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        cubic_eval(xa, ya, &self.c, x, acc)
    }

    fn eval_deriv(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        cubic_eval_deriv(xa, ya, &self.c, x, acc)
    }

    fn eval_deriv2(
        &self,
        xa: &[f64],
        ya: &[f64],
        x: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        cubic_eval_deriv2(xa, ya, &self.c, x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[f64],
        ya: &[f64],
        a: f64,
        b: f64,
        acc: &mut Accelerator,
    ) -> Result<f64, Domain1dError> {
        cubic_eval_integ(xa, ya, &self.c, a, b, acc)
    }
}

//=================================================================================================

fn cubic_eval(
    xa: &[f64],
    ya: &[f64],
    c: &[f64],
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

    let delx = x - xlo;
    let (b, c, d) = coeff_calc(c, dx, dy, index);

    debug_assert!(dx > 0.0);
    Ok(ylo + delx * (b + delx * (c + delx * d)))
}

fn cubic_eval_deriv(
    xa: &[f64],
    ya: &[f64],
    c: &[f64],
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

    let delx = x - xlo;
    let (b, c, d) = coeff_calc(c, dx, dy, index);

    debug_assert!(dx > 0.0);
    Ok(b + delx * (2.0 * c + 3.0 * d * delx))
}

fn cubic_eval_deriv2(
    xa: &[f64],
    ya: &[f64],
    c: &[f64],
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

    let delx = x - xlo;
    let (_, c, d) = coeff_calc(c, dx, dy, index);

    debug_assert!(dx > 0.0);
    Ok(2.0 * c + 6.0 * delx * d)
}

fn cubic_eval_integ(
    xa: &[f64],
    ya: &[f64],
    c: &[f64],
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
        let dy = yhi - ylo;

        // If two x points are the same
        if dx == 0.0 {
            continue;
        }

        let (bi, ci, di) = coeff_calc(c, dx, dy, i);

        if (i == index_a) | (i == index_b) {
            let x1 = if i == index_a { a } else { xlo };
            let x2 = if i == index_b { b } else { xhi };
            result += integ_eval(ylo, bi, ci, di, xlo, x1, x2);
        } else {
            result += dx * (ylo + dx * (0.5 * bi + dx * ((1.0 / 3.0) * ci + 0.25 * di * dx)))
        }
    }
    Ok(result)
}
/// Common coefficient determination.
fn coeff_calc(carray: &[f64], dx: f64, dy: f64, index: usize) -> (f64, f64, f64) {
    let c = carray[index];
    let cplus1 = carray[index + 1];

    let b = (dy / dx) - dx * (cplus1 + 2.0 * c) / 3.0;
    let d = (cplus1 - c) / (3.0 * dx);
    (b, c, d)
}
