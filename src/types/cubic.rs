//! Implementor for Cubic Interpolator.

use std::fmt::Debug;

use ndarray::Array1;
use ndarray_linalg::{Lapack, MatrixLayout, Scalar, SolveTridiagonal, Tridiagonal};
use num::One;

use crate::Accelerator;
use crate::DomainError;
use crate::Interpolation;
use crate::InterpolationError;
use crate::types::utils::integ_eval;
use crate::types::utils::{check_data, check_if_inbounds, diff};

/// Cubic Spline.
///
/// Cubic spline with natural boundary conditions. The resulting curve is piecewise cubic on each
/// interval, with matching first and second derivatives at the supplied data-points. The second
/// derivative is chosen to be zero at the first and last point.
///
/// ## Example
///
/// ```
/// # use rsl_interpolation::Interpolation;
/// # use rsl_interpolation::InterpolationError;
/// # use rsl_interpolation::Cubic;
/// # use rsl_interpolation::Accelerator;
/// #
/// # fn main() -> Result<(), InterpolationError>{
/// let xa = [0.0, 1.0, 2.0];
/// let ya = [0.0, 2.0, 4.0];
/// let interp = Cubic::new(&xa, &ya)?;
/// # Ok(())
/// # }
/// ```
///
/// ## Reference
///
/// Numerical Algorithms with C - Gisela Engeln-Mullges, Frank Uhlig - 1996 -
/// Algorithm 10.1, pg 254
#[allow(dead_code)]
pub struct Cubic<T>
where
    T: num::Float + std::fmt::Debug,
{
    c: Vec<T>,
    g: Vec<T>,
    diag: Vec<T>,
    offdiag: Vec<T>,
}

impl<T> Interpolation<T> for Cubic<T>
where
    T: num::Float + Debug + Scalar + Lapack,
{
    const MIN_SIZE: usize = 3;
    const NAME: &'static str = "cubic";

    fn new(xa: &[T], ya: &[T]) -> Result<Self, crate::InterpolationError>
    where
        Self: Sized,
    {
        check_data(xa, ya, Self::MIN_SIZE)?;

        // Engeln-Mullges G. - Uhlig F.: Algorithm 10.1, pg 254
        let sys_size = xa.len() - 2;

        let h = diff(xa);
        debug_assert!(h.len() == xa.len() - 1);

        let two = T::from(2).unwrap();
        let three = T::from(3).unwrap();

        // Ac=g setup
        let mut g = Vec::<T>::with_capacity(sys_size);
        let mut diag = Vec::<T>::with_capacity(sys_size);
        let mut offdiag = Vec::<T>::with_capacity(sys_size);
        for i in 0..sys_size {
            g.push(if h[i].is_zero() {
                T::zero()
            } else {
                three * (ya[i + 2] - ya[i + 1]) / h[i + 1] - three * (ya[i + 1] - ya[i]) / h[i]
            });
            diag.push(two * (h[i] + h[i + 1]));
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
        let mut c = Vec::<T>::with_capacity(xa.len());
        c.push(T::zero());
        if sys_size.is_one() {
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
        c.push(T::zero());

        // g, diag, and offdiag are only needed for the calculation of c and are not used anywere
        // else from this point, but lets keep them.
        let cubic = Cubic {
            c,
            g,
            diag,
            offdiag,
        };
        Ok(cubic)
    }

    fn eval(&self, xa: &[T], ya: &[T], x: T, acc: &mut Accelerator) -> Result<T, DomainError> {
        cubic_eval(xa, ya, &self.c, x, acc)
    }

    fn eval_deriv(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        cubic_eval_deriv(xa, ya, &self.c, x, acc)
    }

    fn eval_deriv2(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        cubic_eval_deriv2(xa, ya, &self.c, x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[T],
        ya: &[T],
        a: T,
        b: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        cubic_eval_integ(xa, ya, &self.c, a, b, acc)
    }
}

//=================================================================================================

/// Cubic Periodic Spline.
///
/// Cubic Spline with periodic boundary conditions. The resulting curve is piecewise cubic on each
/// interval, with matching first and second derivatives at the supplied data-points. The
/// derivatives at the first and last points are also matched. Note that the last point in the data
/// must have the same y-value as the first point, otherwise the resulting periodic interpolation
/// will have a discontinuity at the boundary.
///
/// ## Example
///
/// ```
/// # use rsl_interpolation::Interpolation;
/// # use rsl_interpolation::InterpolationError;
/// # use rsl_interpolation::CubicPeriodic;
/// # use rsl_interpolation::Accelerator;
/// #
/// # fn main() -> Result<(), InterpolationError>{
/// let xa = [0.0, 1.0, 2.0];
/// let ya = [0.0, 2.0, 4.0];
/// let interp = CubicPeriodic::new(&xa, &ya)?;
/// # Ok(())
/// # }
/// ```
///
/// ## Reference
///
/// Numerical Algorithms with C - Gisela Engeln-Mullges, Frank Uhlig - 1996 -
/// Algorithm 10.2, pg 255
#[allow(dead_code)]
pub struct CubicPeriodic<T>
where
    T: num::Float + std::fmt::Debug,
{
    c: Vec<T>,
    g: Vec<T>,
    diag: Vec<T>,
    offdiag: Vec<T>,
}

impl<T> Interpolation<T> for CubicPeriodic<T>
where
    T: num::Float + Debug + Scalar + Lapack,
{
    const MIN_SIZE: usize = 3;
    const NAME: &'static str = "cubic-periodic";

    fn new(xa: &[T], ya: &[T]) -> Result<Self, crate::InterpolationError>
    where
        Self: Sized,
    {
        check_data(xa, ya, Self::MIN_SIZE)?;

        // Engeln-Mullges G. - Uhlig F.: Algorithm 10.2, pg 255
        let sys_size = xa.len() - 1;

        let h = diff(xa);
        debug_assert!(h.len() == xa.len() - 1);

        let two = T::from(2).unwrap();
        let three = T::from(3).unwrap();

        // Ac=g setup
        let mut c = Vec::<T>::with_capacity(xa.len());
        let mut g = Vec::<T>::with_capacity(sys_size);
        let mut diag = Vec::<T>::with_capacity(sys_size);
        let mut offdiag = Vec::<T>::with_capacity(sys_size);

        if sys_size == 2 {
            let h0 = xa[1] - xa[0];
            let h1 = xa[2] - xa[1];

            let a = two * (h0 + h1);
            let b = h0 + h1;

            g.push(three * ((ya[2] - ya[1]) / h1 - (ya[1] - ya[0]) / h0));
            g.push(three * ((ya[1] - ya[2]) / h0 - (ya[2] - ya[1]) / h1));

            let det = three * (h0 + h1) * (h0 + h1);
            c.push((-b * g[0] + a * g[1]) / det);
            c.push((a * g[0] - b * g[1]) / det);
            c.push(c[0]);
        } else {
            // Same as in Cubic case
            for i in 0..sys_size - 1 {
                g.push(if h[i].is_zero() {
                    T::zero()
                } else {
                    three * (ya[i + 2] - ya[i + 1]) / h[i + 1] - three * (ya[i + 1] - ya[i]) / h[i]
                });
                diag.push(two * (h[i] + h[i + 1]));
                offdiag.push(h[i + 1]);
            }

            // But we must add the last point
            let i = sys_size - 1;
            let hi = xa[i + 1] - xa[i];
            let hiplus1 = xa[1] - xa[0];
            let ydiffi = ya[i + 1] - ya[i];
            let ydiffplus1 = ya[1] - ya[0];
            let gi = if !hi.is_zero() {
                T::one() / hi
            } else {
                T::zero()
            };
            let giplus1 = if !hiplus1.is_zero() {
                T::one() / hiplus1
            } else {
                T::zero()
            };
            offdiag.push(hiplus1);
            diag.push(two * (hiplus1 + hi));
            g.push(three * (ydiffplus1 * giplus1 - ydiffi * gi));
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
            c.push(T::zero());
            if sys_size.is_one() {
                c.push(g[0] / diag[0]);
            } else {
                // This must solve a cyclically tridiagonal matrix, but its not implemented yet :(
                // The corner element is stored at the end of the offdiag vec.
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
            c[0] = c[sys_size];
            panic!(
                "\nNot implemented: Cubic Periodic Splines with more than 3 points require a solver for\
                cyclically tridiagonal matrices, which is currently not implemented by ndarray_linalg.\n"
            )
        }

        // g, diag, and offdiag are only needed for the calculation of c and are not used anywere
        // else from this point, but lets keep them.
        let cubic_periodic = CubicPeriodic {
            c,
            g,
            diag,
            offdiag,
        };
        Ok(cubic_periodic)
    }

    fn eval(&self, xa: &[T], ya: &[T], x: T, acc: &mut Accelerator) -> Result<T, DomainError> {
        cubic_eval(xa, ya, &self.c, x, acc)
    }

    fn eval_deriv(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        cubic_eval_deriv(xa, ya, &self.c, x, acc)
    }

    fn eval_deriv2(
        &self,
        xa: &[T],
        ya: &[T],
        x: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        cubic_eval_deriv2(xa, ya, &self.c, x, acc)
    }

    fn eval_integ(
        &self,
        xa: &[T],
        ya: &[T],
        a: T,
        b: T,
        acc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        cubic_eval_integ(xa, ya, &self.c, a, b, acc)
    }
}

//=================================================================================================

fn cubic_eval<T>(xa: &[T], ya: &[T], c: &[T], x: T, acc: &mut Accelerator) -> Result<T, DomainError>
where
    T: num::Float + Debug + Scalar + Lapack,
{
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

    debug_assert!(dx >= T::zero());
    Ok(ylo + delx * (b + delx * (c + delx * d)))
}

fn cubic_eval_deriv<T>(
    xa: &[T],
    ya: &[T],
    c: &[T],
    x: T,
    acc: &mut Accelerator,
) -> Result<T, DomainError>
where
    T: num::Float + Debug + Scalar + Lapack,
{
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

    let two = T::from(2).unwrap();
    let three = T::from(3).unwrap();

    debug_assert!(dx >= T::zero());
    Ok(b + delx * (two * c + three * d * delx))
}

fn cubic_eval_deriv2<T>(
    xa: &[T],
    ya: &[T],
    c: &[T],
    x: T,
    acc: &mut Accelerator,
) -> Result<T, DomainError>
where
    T: num::Float + Debug + Scalar + Lapack,
{
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

    let two = T::from(2).unwrap();
    let six = T::from(6).unwrap();

    debug_assert!(dx >= T::zero());
    Ok(two * c + six * delx * d)
}

fn cubic_eval_integ<T>(
    xa: &[T],
    ya: &[T],
    c: &[T],
    a: T,
    b: T,
    acc: &mut Accelerator,
) -> Result<T, DomainError>
where
    T: num::Float + Debug + Scalar + Lapack,
{
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
        let ylo = ya[i];
        let yhi = ya[i + 1];

        let dx = xhi - xlo;
        let dy = yhi - ylo;

        // If two x points are the same
        if dx.is_zero() {
            continue;
        }

        let (bi, ci, di) = coeff_calc(c, dx, dy, i);

        if (i == index_a) | (i == index_b) {
            let x1 = if i == index_a { a } else { xlo };
            let x2 = if i == index_b { b } else { xhi };
            result += integ_eval(ylo, bi, ci, di, xlo, x1, x2);
        } else {
            result += dx * (ylo + dx * (half * bi + dx * (third * ci + quarter * di * dx)))
        }
    }
    Ok(result)
}
/// Function for common coefficient determination.
fn coeff_calc<T>(carray: &[T], dx: T, dy: T, index: usize) -> (T, T, T)
where
    T: num::Float + Debug,
{
    let two = T::from(2).unwrap();
    let three = T::from(3).unwrap();

    let c = carray[index];
    let cplus1 = carray[index + 1];

    let b = (dy / dx) - dx * (cplus1 + two * c) / three;
    let d = (cplus1 - c) / (three * dx);
    (b, c, d)
}
