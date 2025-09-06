use ndarray_linalg::Lapack;

use crate::Accelerator;
use crate::Cubic;
use crate::DomainError;
use crate::Interp2dType;
use crate::InterpType;
use crate::Interpolation;
use crate::Interpolation2d;
use crate::InterpolationError;
use crate::interp2d::{acc_indeces, partials, xy_grid_indeces, z_grid_indeces};
use crate::types::utils::check_if_inbounds;
use crate::types::utils::check2d_data;
use crate::z_idx;

const MIN_SIZE: usize = 4;

/// Bicubic Interpolation
#[doc(alias = "gsl_interp2d_bicubic")]
pub struct Bicubic;

impl<T> Interp2dType<T> for Bicubic
where
    T: crate::Num + Lapack,
{
    type Interpolator2d = BicubicInterp<T>;

    /// Constructs a Bicubic Interpolator.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::Interp2dType;
    /// # use rsl_interpolation::Interpolation2d;
    /// # use rsl_interpolation::InterpolationError;
    /// # use rsl_interpolation::Bicubic;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0, 3.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0];
    /// // z = x + y, in column-major order
    /// let za = [
    ///     0.0, 1.0, 2.0, 3.0,
    ///     2.0, 3.0, 4.0, 5.0,
    ///     4.0, 5.0, 6.0, 7.0,
    ///     6.0, 7.0, 8.0, 9.0,
    /// ];
    ///
    /// let interp = Bicubic.build(&xa, &ya, &za)?;
    /// # Ok(())
    /// # }
    /// ```
    fn build(&self, xa: &[T], ya: &[T], za: &[T]) -> Result<BicubicInterp<T>, InterpolationError> {
        check2d_data(xa, ya, za, MIN_SIZE)?;

        let xsize = xa.len();
        let ysize = ya.len();

        // NaN-fill the vecs since their elements are not added in a linear order. NaN ensures that
        // utlimately all coefficients are calculated correctly
        let mut zx: Vec<T> = vec![T::nan(); xsize * ysize];
        let mut zy: Vec<T> = vec![T::nan(); xsize * ysize];
        let mut zxy: Vec<T> = vec![T::nan(); xsize * ysize];

        let mut acc = Accelerator::new();
        let mut x: Vec<T> = vec![T::nan(); xsize];
        let mut y: Vec<T> = vec![T::nan(); xsize];

        #[allow(clippy::needless_range_loop)] // Much cleaner this way
        for j in 0..ysize {
            for i in 0..xsize {
                x[i] = xa[i];
                y[i] = za[z_idx(i, j, xsize, ysize)?];
            }
            let interp = Cubic.build(&x, &y)?;
            for i in 0..xsize {
                let index = z_idx(i, j, xsize, ysize)?;
                zx[index] = interp.eval_deriv(&x, &y, xa[i], &mut acc)?;
            }
        }
        acc.reset(); // Is this necessary?

        let mut x: Vec<T> = vec![T::nan(); ysize];
        let mut y: Vec<T> = vec![T::nan(); ysize];

        #[allow(clippy::needless_range_loop)] // Much cleaner this way
        for i in 0..xsize {
            for j in 0..ysize {
                x[j] = ya[j];
                y[j] = za[z_idx(i, j, xsize, ysize)?];
            }
            let interp = Cubic.build(&x, &y)?;
            for j in 0..ysize {
                let index = z_idx(i, j, xsize, ysize)?;
                zy[index] = interp.eval_deriv(&x, &y, ya[j], &mut acc)?;
            }
        }
        acc.reset();

        let mut x: Vec<T> = vec![T::nan(); xsize];
        let mut y: Vec<T> = vec![T::nan(); xsize];

        #[allow(clippy::needless_range_loop)] // Much cleaner this way
        for j in 0..ysize {
            for i in 0..xsize {
                x[i] = xa[i];
                y[i] = zy[z_idx(i, j, xsize, ysize)?];
            }
            let interp = Cubic.build(&x, &y)?;
            for i in 0..xsize {
                let index = z_idx(i, j, xsize, ysize)?;
                zxy[index] = interp.eval_deriv(&x, &y, xa[i], &mut acc)?;
            }
        }

        let state = BicubicInterp {
            zx,
            zy,
            zxy,
            xsize,
            ysize,
        };

        Ok(state)
    }

    fn name(&self) -> &str {
        "Bicubic"
    }

    fn min_size(&self) -> usize {
        MIN_SIZE
    }
}

// ===============================================================================================

/// Bicubic Interpolator.
///
/// Provides all the evaluation methods.
///
/// Should be constructed through the [`Bicubic`] type.
pub struct BicubicInterp<T>
where
    T: crate::Num + Lapack,
{
    pub(crate) zx: Vec<T>,
    pub(crate) zy: Vec<T>,
    pub(crate) zxy: Vec<T>,
    pub(crate) xsize: usize,
    pub(crate) ysize: usize,
}

impl<T> Interpolation2d<T> for BicubicInterp<T>
where
    T: crate::Num + Lapack,
{
    fn eval_extrap(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        let (xi, yi) = acc_indeces(xa, ya, x, y, xacc, yacc);
        let (xlo, xhi, ylo, yhi) = xy_grid_indeces(xa, ya, xi, yi);
        let (zminmin, zminmax, zmaxmin, zmaxmax) = z_grid_indeces(za, xa.len(), ya.len(), xi, yi)?;
        let (dx, dy) = partials(xlo, xhi, ylo, yhi);

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = self.zxminmaxxing(xi, yi, dt)?;
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = self.zyminmaxxing(xi, yi, du)?;
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = self.zxyminmaxxing(xi, yi, dt, du)?;

        let t2 = t * t;
        let (t0, t1, t2, t3) = (T::one(), t, t2, t * t2);

        let u2 = u * u;
        let (u0, u1, u2, u3) = (T::one(), u, u2, u * u2);

        let two = T::from(2).unwrap();
        let three = T::from(3).unwrap();
        let four = T::from(4).unwrap();
        let six = T::from(6).unwrap();
        let nine = T::from(9).unwrap();

        let mut z = T::zero(); // Result

        let v = zminmin;
        z += v * t0 * u0;
        let v = zyminmin;
        z += v * t0 * u1;
        let v = -three * zminmin + three * zminmax - two * zyminmin - zyminmax;
        z += v * t0 * u2;
        let v = two * zminmin - two * zminmax + zyminmin + zyminmax;
        z += v * t0 * u3;
        let v = zxminmin;
        z += v * t1 * u0;
        let v = zxyminmin;
        z += v * t1 * u1;
        let v = -three * zxminmin + three * zxminmax - two * zxyminmin - zxyminmax;
        z += v * t1 * u2;
        let v = two * zxminmin - two * zxminmax + zxyminmin + zxyminmax;
        z += v * t1 * u3;
        let v = -three * zminmin + three * zmaxmin - two * zxminmin - zxmaxmin;
        z += v * t2 * u0;
        let v = -three * zyminmin + three * zymaxmin - two * zxyminmin - zxymaxmin;
        z += v * t2 * u1;
        #[rustfmt::skip]
        let v = nine * zminmin - nine * zmaxmin + nine * zmaxmax - nine * zminmax
            + six * zxminmin + three * zxmaxmin - three * zxmaxmax - six * zxminmax
            + six * zyminmin - six * zymaxmin - three * zymaxmax + three * zyminmax
            + four * zxyminmin + two * zxymaxmin + zxymaxmax + two * zxyminmax;
        z += v * t2 * u2;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - four * zxminmin - two * zxmaxmin + two * zxmaxmax + four * zxminmax
            - three * zyminmin + three * zymaxmin + three * zymaxmax - three * zyminmax
            - two * zxyminmin - zxymaxmin - zxymaxmax - two * zxyminmax;
        z += v * t2 * u3;
        let v = two * zminmin - two * zmaxmin + zxminmin + zxmaxmin;
        z += v * t3 * u0;
        let v = two * zyminmin - two * zymaxmin + zxyminmin + zxymaxmin;
        z += v * t3 * u1;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - three * zxminmin - three * zxmaxmin + three * zxmaxmax + three * zxminmax
            - four * zyminmin + four * zymaxmin + two * zymaxmax - two * zyminmax
            - two * zxyminmin - two * zxymaxmin - zxymaxmax - zxyminmax;
        z += v * t3 * u2;
        #[rustfmt::skip]
        let v = four * zminmin - four * zmaxmin + four * zmaxmax - four * zminmax
            + two * zxminmin + two * zxmaxmin - two * zxmaxmax - two * zxminmax
            + two * zyminmin - two * zymaxmin - two * zymaxmax + two * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        z += v * t3 * u3;

        Ok(z)
    }

    fn eval_deriv_x(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        check_if_inbounds(xa, x)?;
        check_if_inbounds(ya, y)?;
        let (xi, yi) = acc_indeces(xa, ya, x, y, xacc, yacc);
        let (xlo, xhi, ylo, yhi) = xy_grid_indeces(xa, ya, xi, yi);
        let (zminmin, zminmax, zmaxmin, zmaxmax) = z_grid_indeces(za, xa.len(), ya.len(), xi, yi)?;
        let (dx, dy) = partials(xlo, xhi, ylo, yhi);

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = self.zxminmaxxing(xi, yi, dt)?;
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = self.zyminmaxxing(xi, yi, du)?;
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = self.zxyminmaxxing(xi, yi, dt, du)?;

        let (t0, t1, t2) = (T::one(), t, t * t);

        let u2 = u * u;
        let (u0, u1, u2, u3) = (T::one(), u, u2, u * u2);

        let two = T::from(2).unwrap();
        let three = T::from(3).unwrap();
        let four = T::from(4).unwrap();
        let six = T::from(6).unwrap();
        let nine = T::from(9).unwrap();

        let mut d = T::zero(); // Result

        let v = zxminmin;
        d += v * t0 * u0;
        let v = zxyminmin;
        d += v * t0 * u1;
        let v = -three * zxminmin + three * zxminmax - two * zxyminmin - zxyminmax;
        d += v * t0 * u2;
        let v = two * zxminmin - two * zxminmax + zxyminmin + zxyminmax;
        d += v * t0 * u3;
        let v = -three * zminmin + three * zmaxmin - two * zxminmin - zxmaxmin;
        d += two * v * t1 * u0;
        let v = -three * zyminmin + three * zymaxmin - two * zxyminmin - zxymaxmin;
        d += two * v * t1 * u1;
        #[rustfmt::skip]
        let v = nine * zminmin - nine * zmaxmin + nine * zmaxmax - nine * zminmax
            + six * zxminmin + three * zxmaxmin - three * zxmaxmax - six * zxminmax
            + six * zyminmin - six * zymaxmin - three * zymaxmax + three * zyminmax
            + four * zxyminmin + two * zxymaxmin + zxymaxmax + two * zxyminmax;
        d += two * v * t1 * u2;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - four * zxminmin - two * zxmaxmin + two * zxmaxmax + four * zxminmax
            - three * zyminmin + three * zymaxmin + three * zymaxmax - three * zyminmax
            - two * zxyminmin - zxymaxmin - zxymaxmax - two * zxyminmax;
        d += two * v * t1 * u3;
        let v = two * zminmin - two * zmaxmin + zxminmin + zxmaxmin;
        d += three * v * t2 * u0;
        let v = two * zyminmin - two * zymaxmin + zxyminmin + zxymaxmin;
        d += three * v * t2 * u1;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - three * zxminmin - three * zxmaxmin + three * zxmaxmax + three * zxminmax
            - four * zyminmin + four * zymaxmin + two * zymaxmax - two * zyminmax
            - two * zxyminmin - two * zxymaxmin - zxymaxmax - zxyminmax;
        d += three * v * t2 * u2;
        #[rustfmt::skip]
        let v = four * zminmin - four * zmaxmin + four * zmaxmax - four * zminmax
            + two * zxminmin + two * zxmaxmin - two * zxmaxmax - two * zxminmax
            + two * zyminmin - two * zymaxmin - two * zymaxmax + two * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        d += three * v * t2 * u3;
        d *= dt;

        Ok(d)
    }

    fn eval_deriv_y(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        check_if_inbounds(xa, x)?;
        check_if_inbounds(ya, y)?;
        let (xi, yi) = acc_indeces(xa, ya, x, y, xacc, yacc);
        let (xlo, xhi, ylo, yhi) = xy_grid_indeces(xa, ya, xi, yi);
        let (zminmin, zminmax, zmaxmin, zmaxmax) = z_grid_indeces(za, xa.len(), ya.len(), xi, yi)?;
        let (dx, dy) = partials(xlo, xhi, ylo, yhi);

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = self.zxminmaxxing(xi, yi, dt)?;
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = self.zyminmaxxing(xi, yi, du)?;
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = self.zxyminmaxxing(xi, yi, dt, du)?;

        let t2 = t * t;
        let (t0, t1, t2, t3) = (T::one(), t, t2, t * t2);

        let (u0, u1, u2) = (T::one(), u, u * u);

        let two = T::from(2).unwrap();
        let three = T::from(3).unwrap();
        let four = T::from(4).unwrap();
        let six = T::from(6).unwrap();
        let nine = T::from(9).unwrap();

        let mut d = T::zero(); // Result

        let v = zyminmin;
        d += v * t0 * u0;
        let v = -three * zminmin + three * zminmax - two * zyminmin - zyminmax;
        d += two * v * t0 * u1;
        let v = two * zminmin - two * zminmax + zyminmin + zyminmax;
        d += three * v * t0 * u2;
        let v = zxyminmin;
        d += v * t1 * u0;
        let v = -three * zxminmin + three * zxminmax - two * zxyminmin - zxyminmax;
        d += two * v * t1 * u1;
        let v = two * zxminmin - two * zxminmax + zxyminmin + zxyminmax;
        d += three * v * t1 * u2;
        let v = -three * zyminmin + three * zymaxmin - two * zxyminmin - zxymaxmin;
        d += v * t2 * u0;
        #[rustfmt::skip]
        let v = nine * zminmin - nine * zmaxmin + nine * zmaxmax - nine * zminmax
            + six * zxminmin + three * zxmaxmin - three * zxmaxmax - six * zxminmax
            + six * zyminmin - six * zymaxmin - three * zymaxmax + three * zyminmax
            + four * zxyminmin + two * zxymaxmin + zxymaxmax + two * zxyminmax;
        d += two * v * t2 * u1;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - four * zxminmin - two * zxmaxmin + two * zxmaxmax + four * zxminmax
            - three * zyminmin + three * zymaxmin + three * zymaxmax - three * zyminmax
            - two * zxyminmin - zxymaxmin - zxymaxmax - two * zxyminmax;
        d += three * v * t2 * u2;
        let v = two * zyminmin - two * zymaxmin + zxyminmin + zxymaxmin;
        d += v * t3 * u0;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - three * zxminmin - three * zxmaxmin + three * zxmaxmax + three * zxminmax
            - four * zyminmin + four * zymaxmin + two * zymaxmax - two * zyminmax
            - two * zxyminmin - two * zxymaxmin - zxymaxmax - zxyminmax;
        d += two * v * t3 * u1;
        #[rustfmt::skip]
        let v = four * zminmin - four * zmaxmin + four * zmaxmax - four * zminmax
            + two * zxminmin + two * zxmaxmin - two * zxmaxmax - two * zxminmax
            + two * zyminmin - two * zymaxmin - two * zymaxmax + two * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        d += three * v * t3 * u2;
        d *= du;

        Ok(d)
    }

    fn eval_deriv_xx(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        check_if_inbounds(xa, x)?;
        check_if_inbounds(ya, y)?;
        let (xi, yi) = acc_indeces(xa, ya, x, y, xacc, yacc);
        let (xlo, xhi, ylo, yhi) = xy_grid_indeces(xa, ya, xi, yi);
        let (zminmin, zminmax, zmaxmin, zmaxmax) = z_grid_indeces(za, xa.len(), ya.len(), xi, yi)?;
        let (dx, dy) = partials(xlo, xhi, ylo, yhi);

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = self.zxminmaxxing(xi, yi, dt)?;
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = self.zyminmaxxing(xi, yi, du)?;
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = self.zxyminmaxxing(xi, yi, dt, du)?;

        let (t0, t1) = (T::one(), t);

        let u2 = u * u;
        let (u0, u1, u2, u3) = (T::one(), u, u2, u * u2);

        let two = T::from(2).unwrap();
        let three = T::from(3).unwrap();
        let four = T::from(4).unwrap();
        let six = T::from(6).unwrap();
        let nine = T::from(9).unwrap();

        let mut dd = T::zero(); // Result

        let v = -three * zminmin + three * zmaxmin - two * zxminmin - zxmaxmin;
        dd += two * v * t0 * u0;
        let v = -three * zyminmin + three * zymaxmin - two * zxyminmin - zxymaxmin;
        dd += two * v * t0 * u1;
        #[rustfmt::skip]
        let v = nine * zminmin - nine * zmaxmin + nine * zmaxmax - nine * zminmax
            + six * zxminmin + three * zxmaxmin - three * zxmaxmax - six * zxminmax
            + six * zyminmin - six * zymaxmin - three * zymaxmax + three * zyminmax
            + four * zxyminmin + two * zxymaxmin + zxymaxmax + two * zxyminmax;
        dd += two * v * t0 * u2;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - four * zxminmin - two * zxmaxmin + two * zxmaxmax + four * zxminmax
            - three * zyminmin + three * zymaxmin + three * zymaxmax - three * zyminmax
            - two * zxyminmin - zxymaxmin - zxymaxmax - two * zxyminmax;
        dd += two * v * t0 * u3;
        let v = two * zminmin - two * zmaxmin + zxminmin + zxmaxmin;
        dd += six * v * t1 * u0;
        let v = two * zyminmin - two * zymaxmin + zxyminmin + zxymaxmin;
        dd += six * v * t1 * u1;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - three * zxminmin - three * zxmaxmin + three * zxmaxmax + three * zxminmax
            - four * zyminmin + four * zymaxmin + two * zymaxmax - two * zyminmax
            - two * zxyminmin - two * zxymaxmin - zxymaxmax - zxyminmax;
        dd += six * v * t1 * u2;
        #[rustfmt::skip]
        let v = four * zminmin - four * zmaxmin + four * zmaxmax - four * zminmax
            + two * zxminmin + two * zxmaxmin - two * zxmaxmax - two * zxminmax
            + two * zyminmin - two * zymaxmin - two * zymaxmax + two * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        dd += six * v * t1 * u3;
        dd = dd * dt * dt;

        Ok(dd)
    }

    fn eval_deriv_yy(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        check_if_inbounds(xa, x)?;
        check_if_inbounds(ya, y)?;
        let (xi, yi) = acc_indeces(xa, ya, x, y, xacc, yacc);
        let (xlo, xhi, ylo, yhi) = xy_grid_indeces(xa, ya, xi, yi);
        let (zminmin, zminmax, zmaxmin, zmaxmax) = z_grid_indeces(za, xa.len(), ya.len(), xi, yi)?;
        let (dx, dy) = partials(xlo, xhi, ylo, yhi);

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = self.zxminmaxxing(xi, yi, dt)?;
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = self.zyminmaxxing(xi, yi, du)?;
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = self.zxyminmaxxing(xi, yi, dt, du)?;

        let t2 = t * t;
        let (t0, t1, t2, t3) = (T::one(), t, t2, t * t2);

        let (u0, u1) = (T::one(), u);

        let two = T::from(2).unwrap();
        let three = T::from(3).unwrap();
        let four = T::from(4).unwrap();
        let six = T::from(6).unwrap();
        let nine = T::from(9).unwrap();

        let mut dd = T::zero(); // Result

        let v = -three * zminmin + three * zminmax - two * zyminmin - zyminmax;
        dd += two * v * t0 * u0;
        let v = two * zminmin - two * zminmax + zyminmin + zyminmax;
        dd += six * v * t0 * u1;
        let v = -three * zxminmin + three * zxminmax - two * zxyminmin - zxyminmax;
        dd += two * v * t1 * u0;
        let v = two * zxminmin - two * zxminmax + zxyminmin + zxyminmax;
        dd += six * v * t1 * u1;
        #[rustfmt::skip]
        let v = nine * zminmin - nine * zmaxmin + nine * zmaxmax - nine * zminmax
            + six * zxminmin + three * zxmaxmin - three * zxmaxmax - six * zxminmax
            + six * zyminmin - six * zymaxmin - three * zymaxmax + three * zyminmax
            + four * zxyminmin + two * zxymaxmin + zxymaxmax + two * zxyminmax;
        dd += two * v * t2 * u0;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - four * zxminmin - two * zxmaxmin + two * zxmaxmax + four * zxminmax
            - three * zyminmin + three * zymaxmin + three * zymaxmax - three * zyminmax
            - two * zxyminmin - zxymaxmin - zxymaxmax - two * zxyminmax;
        dd += six * v * t2 * u1;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - three * zxminmin - three * zxmaxmin + three * zxmaxmax + three * zxminmax
            - four * zyminmin + four * zymaxmin + two * zymaxmax - two * zyminmax
            - two * zxyminmin - two * zxymaxmin - zxymaxmax - zxyminmax;
        dd += two * v * t3 * u0;
        #[rustfmt::skip]
        let v = four * zminmin - four * zmaxmin + four * zmaxmax - four * zminmax
            + two * zxminmin + two * zxmaxmin - two * zxmaxmax - two * zxminmax
            + two * zyminmin - two * zymaxmin - two * zymaxmax + two * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        dd += six * v * t3 * u1;
        dd = dd * du * du;

        Ok(dd)
    }

    fn eval_deriv_xy(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        check_if_inbounds(xa, x)?;
        check_if_inbounds(ya, y)?;
        let (xi, yi) = acc_indeces(xa, ya, x, y, xacc, yacc);
        let (xlo, xhi, ylo, yhi) = xy_grid_indeces(xa, ya, xi, yi);
        let (zminmin, zminmax, zmaxmin, zmaxmax) = z_grid_indeces(za, xa.len(), ya.len(), xi, yi)?;
        let (dx, dy) = partials(xlo, xhi, ylo, yhi);

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = self.zxminmaxxing(xi, yi, dt)?;
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = self.zyminmaxxing(xi, yi, du)?;
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = self.zxyminmaxxing(xi, yi, dt, du)?;

        let (t0, t1, t2) = (T::one(), t, t * t);

        let (u0, u1, u2) = (T::one(), u, u * u);

        let two = T::from(2).unwrap();
        let three = T::from(3).unwrap();
        let four = T::from(4).unwrap();
        let six = T::from(6).unwrap();
        let nine = T::from(9).unwrap();

        let mut dd = T::zero(); // Result

        let v = zxyminmin;
        dd += v * t0 * u0;
        let v = -three * zxminmin + three * zxminmax - two * zxyminmin - zxyminmax;
        dd += two * v * t0 * u1;
        let v = two * zxminmin - two * zxminmax + zxyminmin + zxyminmax;
        dd += three * v * t0 * u2;
        let v = -three * zyminmin + three * zymaxmin - two * zxyminmin - zxymaxmin;
        dd += two * v * t1 * u0;
        #[rustfmt::skip]
        let v = nine * zminmin - nine * zmaxmin + nine * zmaxmax - nine * zminmax
            + six * zxminmin + three * zxmaxmin - three * zxmaxmax - six * zxminmax
            + six * zyminmin - six * zymaxmin - three * zymaxmax + three * zyminmax
            + four * zxyminmin + two * zxymaxmin + zxymaxmax + two * zxyminmax;
        dd += four * v * t1 * u1;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - four * zxminmin - two * zxmaxmin + two * zxmaxmax + four * zxminmax
            - three * zyminmin + three * zymaxmin + three * zymaxmax - three * zyminmax
            - two * zxyminmin - zxymaxmin - zxymaxmax - two * zxyminmax;
        dd += six * v * t1 * u2;
        let v = two * zyminmin - two * zymaxmin + zxyminmin + zxymaxmin;
        dd += three * v * t2 * u0;
        #[rustfmt::skip]
        let v = -six * zminmin + six * zmaxmin - six * zmaxmax + six * zminmax
            - three * zxminmin - three * zxmaxmin + three * zxmaxmax + three * zxminmax
            - four * zyminmin + four * zymaxmin + two * zymaxmax - two * zyminmax
            - two * zxyminmin - two * zxymaxmin - zxymaxmax - zxyminmax;
        dd += six * v * t2 * u1;
        #[rustfmt::skip]
        let v = four * zminmin - four * zmaxmin + four * zmaxmax - four * zminmax
            + two * zxminmin + two * zxmaxmin - two * zxmaxmax - two * zxminmax
            + two * zyminmin - two * zymaxmin - two * zymaxmax + two * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        dd += nine * v * t2 * u2;
        dd = dd * dt * du;

        Ok(dd)
    }
}

/// Common calculations
impl<T> BicubicInterp<T>
where
    T: crate::Num + Lapack,
{
    fn zxminmaxxing(&self, xi: usize, yi: usize, dt: T) -> Result<(T, T, T, T), DomainError>
    where
        T: crate::Num + Lapack,
    {
        let zxminmin = self.zx[z_idx(xi, yi, self.xsize, self.ysize)?] / dt;
        let zxminmax = self.zx[z_idx(xi, yi + 1, self.xsize, self.ysize)?] / dt;
        let zxmaxmin = self.zx[z_idx(xi + 1, yi, self.xsize, self.ysize)?] / dt;
        let zxmaxmax = self.zx[z_idx(xi + 1, yi + 1, self.xsize, self.ysize)?] / dt;
        Ok((zxminmin, zxminmax, zxmaxmin, zxmaxmax))
    }

    fn zyminmaxxing(&self, xi: usize, yi: usize, du: T) -> Result<(T, T, T, T), DomainError>
    where
        T: crate::Num + Lapack,
    {
        let zyminmin = self.zy[z_idx(xi, yi, self.xsize, self.ysize)?] / du;
        let zyminmax = self.zy[z_idx(xi, yi + 1, self.xsize, self.ysize)?] / du;
        let zymaxmin = self.zy[z_idx(xi + 1, yi, self.xsize, self.ysize)?] / du;
        let zymaxmax = self.zy[z_idx(xi + 1, yi + 1, self.xsize, self.ysize)?] / du;
        Ok((zyminmin, zyminmax, zymaxmin, zymaxmax))
    }

    fn zxyminmaxxing(&self, xi: usize, yi: usize, dt: T, du: T) -> Result<(T, T, T, T), DomainError>
    where
        T: crate::Num + Lapack,
    {
        let prod = dt * du;
        let zxyminmin = self.zxy[z_idx(xi, yi, self.xsize, self.ysize)?] / prod;
        let zxyminmax = self.zxy[z_idx(xi, yi + 1, self.xsize, self.ysize)?] / prod;
        let zxymaxmin = self.zxy[z_idx(xi + 1, yi, self.xsize, self.ysize)?] / prod;
        let zxymaxmax = self.zxy[z_idx(xi + 1, yi + 1, self.xsize, self.ysize)?] / prod;
        Ok((zxyminmin, zxyminmax, zxymaxmin, zxymaxmax))
    }
}

/// Common calculation
fn tu_cubic_values<T>(x: T, y: T, xlo: T, ylo: T, dx: T, dy: T) -> (T, T, T, T)
where
    T: crate::Num + Lapack,
{
    let t = (x - xlo) / dx;
    let u = (y - ylo) / dy;
    let dt = dx.recip();
    let du = dy.recip();
    (t, u, dt, du)
}
