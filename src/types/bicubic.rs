//! Definition of `Bicubic` Interpolator.

use crate::{
    Accelerator, Accelerator2d, CubicInterpolator, Domain2dError, Interpolation, Interpolation2d,
    InterpolationError,
};
use crate::{check_if_inbounds2d, check2d_data, z_idx};

/// Bicubic Interpolator.
#[doc(alias = "gsl_interp2d_bicubic")]
#[derive(Debug, Clone)]
pub struct BicubicInterpolator {
    pub(crate) zx: Box<[f64]>,
    pub(crate) zy: Box<[f64]>,
    pub(crate) zxy: Box<[f64]>,
}

impl BicubicInterpolator {
    /// The minimum required number of points.
    #[doc(alias = "gsl_interp2d_min_size")]
    const MIN_SIZE: usize = 4;

    /// Constructs a Bicubic Interpolator.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
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
    /// let interp = BicubicInterpolator::build(&xa, &ya, &za)?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// - [`InterpolationError::UnsortedDataset`]: `xa` or `ya` are not monotonically increasing.
    /// - [`InterpolationError::NotEnoughPoints`]: length of `xa` or `ya` is less that 4.
    /// - [`InterpolationError::BLASTridiagError`]: Error when solving the tridiagonal system.
    /// - [`InterpolationError::ZGridMismatch`]: `xa.len()*ya.len() != za.len()`.
    #[expect(clippy::needless_range_loop, reason = "Much cleaner this way")]
    #[expect(clippy::missing_panics_doc, reason = "`xa` and `ya` are checked")]
    #[expect(clippy::unwrap_in_result, reason = "`xa` and `ya` are checked")]
    #[doc(alias = "gsl_interp2d_init")]
    pub fn build(xa: &[f64], ya: &[f64], za: &[f64]) -> Result<Self, InterpolationError> {
        check2d_data(xa, ya, za, Self::MIN_SIZE)?;

        let xsize = xa.len();
        let ysize = ya.len();

        // NaN-fill the vecs since their elements are not added in a linear order. NaN ensures that
        // ultimately all coefficients are calculated correctly
        let mut zx = vec![f64::NAN; xsize * ysize];
        let mut zy = vec![f64::NAN; xsize * ysize];
        let mut zxy = vec![f64::NAN; xsize * ysize];

        let mut acc = Accelerator::new();
        let mut x = vec![f64::NAN; xsize];
        let mut y = vec![f64::NAN; xsize];

        for j in 0..ysize {
            for i in 0..xsize {
                x[i] = xa[i];
                y[i] = za[z_idx(i, j, xsize, ysize)];
            }
            let interp = CubicInterpolator::build(&x, &y).expect("checked");
            for i in 0..xsize {
                let index = z_idx(i, j, xsize, ysize);
                zx[index] = interp.eval_deriv(&x, &y, xa[i], &mut acc).expect("checked");
            }
        }

        acc.reset();
        let mut x = vec![f64::NAN; ysize];
        let mut y = vec![f64::NAN; ysize];

        for i in 0..xsize {
            for j in 0..ysize {
                x[j] = ya[j];
                y[j] = za[z_idx(i, j, xsize, ysize)];
            }
            let interp = CubicInterpolator::build(&x, &y).expect("checked");
            for j in 0..ysize {
                let index = z_idx(i, j, xsize, ysize);
                zy[index] = interp.eval_deriv(&x, &y, ya[j], &mut acc).expect("checked");
            }
        }

        acc.reset();
        let mut x = vec![f64::NAN; xsize];
        let mut y = vec![f64::NAN; xsize];

        for j in 0..ysize {
            for i in 0..xsize {
                x[i] = xa[i];
                y[i] = zy[z_idx(i, j, xsize, ysize)];
            }
            let interp = CubicInterpolator::build(&x, &y).expect("checked");
            for i in 0..xsize {
                let index = z_idx(i, j, xsize, ysize);
                zxy[index] = interp.eval_deriv(&x, &y, xa[i], &mut acc).expect("checked");
            }
        }

        Ok(BicubicInterpolator {
            zx: zx.into_boxed_slice(),
            zy: zy.into_boxed_slice(),
            zxy: zxy.into_boxed_slice(),
        })
    }
}

// ===============================================================================================

impl Interpolation2d for BicubicInterpolator {
    fn eval_extrap(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> f64 {
        let is_uptodate = acc.is_uptodate(xa, ya, x, y);
        if !is_uptodate {
            acc.update_step1(xa, ya, za);
        }

        let (xlo, _, ylo, _) = acc.get_xy_grid_values();
        let (zminmin, zminmax, zmaxmin, zmaxmax) = acc.get_z_grid_values();
        let (dx, dy) = acc.get_partials();

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        if !is_uptodate {
            acc.update_step2(xa, ya, &self.zx, &self.zy, &self.zxy, dt, du);
        };

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = acc.get_zxminmaxxing();
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = acc.get_zyminmaxxing();
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = acc.get_zxyminmaxxing();

        let t2 = t * t;
        let (t0, t1, t2, t3) = (1.0, t, t2, t * t2);

        let u2 = u * u;
        let (u0, u1, u2, u3) = (1.0, u, u2, u * u2);

        let mut z = 0.0; // Result

        let v = zminmin;
        z += v * t0 * u0;
        let v = zyminmin;
        z += v * t0 * u1;
        let v = -3.0 * zminmin + 3.0 * zminmax - 2.0 * zyminmin - zyminmax;
        z += v * t0 * u2;
        let v = 2.0 * zminmin - 2.0 * zminmax + zyminmin + zyminmax;
        z += v * t0 * u3;
        let v = zxminmin;
        z += v * t1 * u0;
        let v = zxyminmin;
        z += v * t1 * u1;
        let v = -3.0 * zxminmin + 3.0 * zxminmax - 2.0 * zxyminmin - zxyminmax;
        z += v * t1 * u2;
        let v = 2.0 * zxminmin - 2.0 * zxminmax + zxyminmin + zxyminmax;
        z += v * t1 * u3;
        let v = -3.0 * zminmin + 3.0 * zmaxmin - 2.0 * zxminmin - zxmaxmin;
        z += v * t2 * u0;
        let v = -3.0 * zyminmin + 3.0 * zymaxmin - 2.0 * zxyminmin - zxymaxmin;
        z += v * t2 * u1;
        #[rustfmt::skip]
        let v = 9.0 * zminmin - 9.0 * zmaxmin + 9.0 * zmaxmax - 9.0 * zminmax
            + 6.0 * zxminmin + 3.0 * zxmaxmin - 3.0 * zxmaxmax - 6.0 * zxminmax
            + 6.0 * zyminmin - 6.0 * zymaxmin - 3.0 * zymaxmax + 3.0 * zyminmax
            + 4.0 * zxyminmin + 2.0 * zxymaxmin + zxymaxmax + 2.0 * zxyminmax;
        z += v * t2 * u2;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 4.0 * zxminmin - 2.0 * zxmaxmin + 2.0 * zxmaxmax + 4.0 * zxminmax
            - 3.0 * zyminmin + 3.0 * zymaxmin + 3.0 * zymaxmax - 3.0 * zyminmax
            - 2.0 * zxyminmin - zxymaxmin - zxymaxmax - 2.0 * zxyminmax;
        z += v * t2 * u3;
        let v = 2.0 * zminmin - 2.0 * zmaxmin + zxminmin + zxmaxmin;
        z += v * t3 * u0;
        let v = 2.0 * zyminmin - 2.0 * zymaxmin + zxyminmin + zxymaxmin;
        z += v * t3 * u1;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 3.0 * zxminmin - 3.0 * zxmaxmin + 3.0 * zxmaxmax + 3.0 * zxminmax
            - 4.0 * zyminmin + 4.0 * zymaxmin + 2.0 * zymaxmax - 2.0 * zyminmax
            - 2.0 * zxyminmin - 2.0 * zxymaxmin - zxymaxmax - zxyminmax;
        z += v * t3 * u2;
        #[rustfmt::skip]
        let v = 4.0 * zminmin - 4.0 * zmaxmin + 4.0 * zmaxmax - 4.0 * zminmax
            + 2.0 * zxminmin + 2.0 * zxmaxmin - 2.0 * zxmaxmax - 2.0 * zxminmax
            + 2.0 * zyminmin - 2.0 * zymaxmin - 2.0 * zymaxmax + 2.0 * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        z += v * t3 * u3;

        z
    }

    fn eval_deriv_x(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        check_if_inbounds2d(xa, ya, x, y)?;
        let is_uptodate = acc.is_uptodate(xa, ya, x, y);
        if !is_uptodate {
            acc.update_step1(xa, ya, za);
        }

        let (xlo, _, ylo, _) = acc.get_xy_grid_values();
        let (zminmin, zminmax, zmaxmin, zmaxmax) = acc.get_z_grid_values();
        let (dx, dy) = acc.get_partials();

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        if !is_uptodate {
            acc.update_step2(xa, ya, &self.zx, &self.zy, &self.zxy, dt, du);
        };

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = acc.get_zxminmaxxing();
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = acc.get_zyminmaxxing();
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = acc.get_zxyminmaxxing();

        let (t0, t1, t2) = (1.0, t, t * t);

        let u2 = u * u;
        let (u0, u1, u2, u3) = (1.0, u, u2, u * u2);

        let mut d = 0.0; // Result

        let v = zxminmin;
        d += v * t0 * u0;
        let v = zxyminmin;
        d += v * t0 * u1;
        let v = -3.0 * zxminmin + 3.0 * zxminmax - 2.0 * zxyminmin - zxyminmax;
        d += v * t0 * u2;
        let v = 2.0 * zxminmin - 2.0 * zxminmax + zxyminmin + zxyminmax;
        d += v * t0 * u3;
        let v = -3.0 * zminmin + 3.0 * zmaxmin - 2.0 * zxminmin - zxmaxmin;
        d += 2.0 * v * t1 * u0;
        let v = -3.0 * zyminmin + 3.0 * zymaxmin - 2.0 * zxyminmin - zxymaxmin;
        d += 2.0 * v * t1 * u1;
        #[rustfmt::skip]
        let v = 9.0 * zminmin - 9.0 * zmaxmin + 9.0 * zmaxmax - 9.0 * zminmax
            + 6.0 * zxminmin + 3.0 * zxmaxmin - 3.0 * zxmaxmax - 6.0 * zxminmax
            + 6.0 * zyminmin - 6.0 * zymaxmin - 3.0 * zymaxmax + 3.0 * zyminmax
            + 4.0 * zxyminmin + 2.0 * zxymaxmin + zxymaxmax + 2.0 * zxyminmax;
        d += 2.0 * v * t1 * u2;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 4.0 * zxminmin - 2.0 * zxmaxmin + 2.0 * zxmaxmax + 4.0 * zxminmax
            - 3.0 * zyminmin + 3.0 * zymaxmin + 3.0 * zymaxmax - 3.0 * zyminmax
            - 2.0 * zxyminmin - zxymaxmin - zxymaxmax - 2.0 * zxyminmax;
        d += 2.0 * v * t1 * u3;
        let v = 2.0 * zminmin - 2.0 * zmaxmin + zxminmin + zxmaxmin;
        d += 3.0 * v * t2 * u0;
        let v = 2.0 * zyminmin - 2.0 * zymaxmin + zxyminmin + zxymaxmin;
        d += 3.0 * v * t2 * u1;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 3.0 * zxminmin - 3.0 * zxmaxmin + 3.0 * zxmaxmax + 3.0 * zxminmax
            - 4.0 * zyminmin + 4.0 * zymaxmin + 2.0 * zymaxmax - 2.0 * zyminmax
            - 2.0 * zxyminmin - 2.0 * zxymaxmin - zxymaxmax - zxyminmax;
        d += 3.0 * v * t2 * u2;
        #[rustfmt::skip]
        let v = 4.0 * zminmin - 4.0 * zmaxmin + 4.0 * zmaxmax - 4.0 * zminmax
            + 2.0 * zxminmin + 2.0 * zxmaxmin - 2.0 * zxmaxmax - 2.0 * zxminmax
            + 2.0 * zyminmin - 2.0 * zymaxmin - 2.0 * zymaxmax + 2.0 * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        d += 3.0 * v * t2 * u3;
        d *= dt;

        Ok(d)
    }

    fn eval_deriv_y(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        check_if_inbounds2d(xa, ya, x, y)?;
        let is_uptodate = acc.is_uptodate(xa, ya, x, y);
        if !is_uptodate {
            acc.update_step1(xa, ya, za);
        }

        let (xlo, _, ylo, _) = acc.get_xy_grid_values();
        let (zminmin, zminmax, zmaxmin, zmaxmax) = acc.get_z_grid_values();
        let (dx, dy) = acc.get_partials();

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        if !is_uptodate {
            acc.update_step2(xa, ya, &self.zx, &self.zy, &self.zxy, dt, du);
        };

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = acc.get_zxminmaxxing();
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = acc.get_zyminmaxxing();
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = acc.get_zxyminmaxxing();

        let t2 = t * t;
        let (t0, t1, t2, t3) = (1.0, t, t2, t * t2);

        let (u0, u1, u2) = (1.0, u, u * u);

        let mut d = 0.0; // Result

        let v = zyminmin;
        d += v * t0 * u0;
        let v = -3.0 * zminmin + 3.0 * zminmax - 2.0 * zyminmin - zyminmax;
        d += 2.0 * v * t0 * u1;
        let v = 2.0 * zminmin - 2.0 * zminmax + zyminmin + zyminmax;
        d += 3.0 * v * t0 * u2;
        let v = zxyminmin;
        d += v * t1 * u0;
        let v = -3.0 * zxminmin + 3.0 * zxminmax - 2.0 * zxyminmin - zxyminmax;
        d += 2.0 * v * t1 * u1;
        let v = 2.0 * zxminmin - 2.0 * zxminmax + zxyminmin + zxyminmax;
        d += 3.0 * v * t1 * u2;
        let v = -3.0 * zyminmin + 3.0 * zymaxmin - 2.0 * zxyminmin - zxymaxmin;
        d += v * t2 * u0;
        #[rustfmt::skip]
        let v = 9.0 * zminmin - 9.0 * zmaxmin + 9.0 * zmaxmax - 9.0 * zminmax
            + 6.0 * zxminmin + 3.0 * zxmaxmin - 3.0 * zxmaxmax - 6.0 * zxminmax
            + 6.0 * zyminmin - 6.0 * zymaxmin - 3.0 * zymaxmax + 3.0 * zyminmax
            + 4.0 * zxyminmin + 2.0 * zxymaxmin + zxymaxmax + 2.0 * zxyminmax;
        d += 2.0 * v * t2 * u1;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 4.0 * zxminmin - 2.0 * zxmaxmin + 2.0 * zxmaxmax + 4.0 * zxminmax
            - 3.0 * zyminmin + 3.0 * zymaxmin + 3.0 * zymaxmax - 3.0 * zyminmax
            - 2.0 * zxyminmin - zxymaxmin - zxymaxmax - 2.0 * zxyminmax;
        d += 3.0 * v * t2 * u2;
        let v = 2.0 * zyminmin - 2.0 * zymaxmin + zxyminmin + zxymaxmin;
        d += v * t3 * u0;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 3.0 * zxminmin - 3.0 * zxmaxmin + 3.0 * zxmaxmax + 3.0 * zxminmax
            - 4.0 * zyminmin + 4.0 * zymaxmin + 2.0 * zymaxmax - 2.0 * zyminmax
            - 2.0 * zxyminmin - 2.0 * zxymaxmin - zxymaxmax - zxyminmax;
        d += 2.0 * v * t3 * u1;
        #[rustfmt::skip]
        let v = 4.0 * zminmin - 4.0 * zmaxmin + 4.0 * zmaxmax - 4.0 * zminmax
            + 2.0 * zxminmin + 2.0 * zxmaxmin - 2.0 * zxmaxmax - 2.0 * zxminmax
            + 2.0 * zyminmin - 2.0 * zymaxmin - 2.0 * zymaxmax + 2.0 * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        d += 3.0 * v * t3 * u2;
        d *= du;

        Ok(d)
    }

    fn eval_deriv_xx(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        check_if_inbounds2d(xa, ya, x, y)?;
        let is_uptodate = acc.is_uptodate(xa, ya, x, y);
        if !is_uptodate {
            acc.update_step1(xa, ya, za);
        }

        let (xlo, _, ylo, _) = acc.get_xy_grid_values();
        let (zminmin, zminmax, zmaxmin, zmaxmax) = acc.get_z_grid_values();
        let (dx, dy) = acc.get_partials();

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        if !is_uptodate {
            acc.update_step2(xa, ya, &self.zx, &self.zy, &self.zxy, dt, du);
        };

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = acc.get_zxminmaxxing();
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = acc.get_zyminmaxxing();
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = acc.get_zxyminmaxxing();

        let (t0, t1) = (1.0, t);

        let u2 = u * u;
        let (u0, u1, u2, u3) = (1.0, u, u2, u * u2);

        let mut dd = 0.0; // Result

        let v = -3.0 * zminmin + 3.0 * zmaxmin - 2.0 * zxminmin - zxmaxmin;
        dd += 2.0 * v * t0 * u0;
        let v = -3.0 * zyminmin + 3.0 * zymaxmin - 2.0 * zxyminmin - zxymaxmin;
        dd += 2.0 * v * t0 * u1;
        #[rustfmt::skip]
        let v = 9.0 * zminmin - 9.0 * zmaxmin + 9.0 * zmaxmax - 9.0 * zminmax
            + 6.0 * zxminmin + 3.0 * zxmaxmin - 3.0 * zxmaxmax - 6.0 * zxminmax
            + 6.0 * zyminmin - 6.0 * zymaxmin - 3.0 * zymaxmax + 3.0 * zyminmax
            + 4.0 * zxyminmin + 2.0 * zxymaxmin + zxymaxmax + 2.0 * zxyminmax;
        dd += 2.0 * v * t0 * u2;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 4.0 * zxminmin - 2.0 * zxmaxmin + 2.0 * zxmaxmax + 4.0 * zxminmax
            - 3.0 * zyminmin + 3.0 * zymaxmin + 3.0 * zymaxmax - 3.0 * zyminmax
            - 2.0 * zxyminmin - zxymaxmin - zxymaxmax - 2.0 * zxyminmax;
        dd += 2.0 * v * t0 * u3;
        let v = 2.0 * zminmin - 2.0 * zmaxmin + zxminmin + zxmaxmin;
        dd += 6.0 * v * t1 * u0;
        let v = 2.0 * zyminmin - 2.0 * zymaxmin + zxyminmin + zxymaxmin;
        dd += 6.0 * v * t1 * u1;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 3.0 * zxminmin - 3.0 * zxmaxmin + 3.0 * zxmaxmax + 3.0 * zxminmax
            - 4.0 * zyminmin + 4.0 * zymaxmin + 2.0 * zymaxmax - 2.0 * zyminmax
            - 2.0 * zxyminmin - 2.0 * zxymaxmin - zxymaxmax - zxyminmax;
        dd += 6.0 * v * t1 * u2;
        #[rustfmt::skip]
        let v = 4.0 * zminmin - 4.0 * zmaxmin + 4.0 * zmaxmax - 4.0 * zminmax
            + 2.0 * zxminmin + 2.0 * zxmaxmin - 2.0 * zxmaxmax - 2.0 * zxminmax
            + 2.0 * zyminmin - 2.0 * zymaxmin - 2.0 * zymaxmax + 2.0 * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        dd += 6.0 * v * t1 * u3;
        dd = dd * dt * dt;

        Ok(dd)
    }

    fn eval_deriv_yy(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        check_if_inbounds2d(xa, ya, x, y)?;
        let is_uptodate = acc.is_uptodate(xa, ya, x, y);
        if !is_uptodate {
            acc.update_step1(xa, ya, za);
        }

        let (xlo, _, ylo, _) = acc.get_xy_grid_values();
        let (zminmin, zminmax, zmaxmin, zmaxmax) = acc.get_z_grid_values();
        let (dx, dy) = acc.get_partials();

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        if !is_uptodate {
            acc.update_step2(xa, ya, &self.zx, &self.zy, &self.zxy, dt, du);
        };

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = acc.get_zxminmaxxing();
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = acc.get_zyminmaxxing();
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = acc.get_zxyminmaxxing();

        let t2 = t * t;
        let (t0, t1, t2, t3) = (1.0, t, t2, t * t2);

        let (u0, u1) = (1.0, u);

        let mut dd = 0.0; // Result

        let v = -3.0 * zminmin + 3.0 * zminmax - 2.0 * zyminmin - zyminmax;
        dd += 2.0 * v * t0 * u0;
        let v = 2.0 * zminmin - 2.0 * zminmax + zyminmin + zyminmax;
        dd += 6.0 * v * t0 * u1;
        let v = -3.0 * zxminmin + 3.0 * zxminmax - 2.0 * zxyminmin - zxyminmax;
        dd += 2.0 * v * t1 * u0;
        let v = 2.0 * zxminmin - 2.0 * zxminmax + zxyminmin + zxyminmax;
        dd += 6.0 * v * t1 * u1;
        #[rustfmt::skip]
        let v = 9.0 * zminmin - 9.0 * zmaxmin + 9.0 * zmaxmax - 9.0 * zminmax
            + 6.0 * zxminmin + 3.0 * zxmaxmin - 3.0 * zxmaxmax - 6.0 * zxminmax
            + 6.0 * zyminmin - 6.0 * zymaxmin - 3.0 * zymaxmax + 3.0 * zyminmax
            + 4.0 * zxyminmin + 2.0 * zxymaxmin + zxymaxmax + 2.0 * zxyminmax;
        dd += 2.0 * v * t2 * u0;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 4.0 * zxminmin - 2.0 * zxmaxmin + 2.0 * zxmaxmax + 4.0 * zxminmax
            - 3.0 * zyminmin + 3.0 * zymaxmin + 3.0 * zymaxmax - 3.0 * zyminmax
            - 2.0 * zxyminmin - zxymaxmin - zxymaxmax - 2.0 * zxyminmax;
        dd += 6.0 * v * t2 * u1;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 3.0 * zxminmin - 3.0 * zxmaxmin + 3.0 * zxmaxmax + 3.0 * zxminmax
            - 4.0 * zyminmin + 4.0 * zymaxmin + 2.0 * zymaxmax - 2.0 * zyminmax
            - 2.0 * zxyminmin - 2.0 * zxymaxmin - zxymaxmax - zxyminmax;
        dd += 2.0 * v * t3 * u0;
        #[rustfmt::skip]
        let v = 4.0 * zminmin - 4.0 * zmaxmin + 4.0 * zmaxmax - 4.0 * zminmax
            + 2.0 * zxminmin + 2.0 * zxmaxmin - 2.0 * zxmaxmax - 2.0 * zxminmax
            + 2.0 * zyminmin - 2.0 * zymaxmin - 2.0 * zymaxmax + 2.0 * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        dd += 6.0 * v * t3 * u1;
        dd = dd * du * du;

        Ok(dd)
    }

    fn eval_deriv_xy(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        check_if_inbounds2d(xa, ya, x, y)?;
        let is_uptodate = acc.is_uptodate(xa, ya, x, y);
        if !is_uptodate {
            acc.update_step1(xa, ya, za);
        }

        let (xlo, _, ylo, _) = acc.get_xy_grid_values();
        let (zminmin, zminmax, zmaxmin, zmaxmax) = acc.get_z_grid_values();
        let (dx, dy) = acc.get_partials();

        let (t, u, dt, du) = tu_cubic_values(x, y, xlo, ylo, dx, dy);

        if !is_uptodate {
            acc.update_step2(xa, ya, &self.zx, &self.zy, &self.zxy, dt, du);
        };

        let (zxminmin, zxminmax, zxmaxmin, zxmaxmax) = acc.get_zxminmaxxing();
        let (zyminmin, zyminmax, zymaxmin, zymaxmax) = acc.get_zyminmaxxing();
        let (zxyminmin, zxyminmax, zxymaxmin, zxymaxmax) = acc.get_zxyminmaxxing();

        let (t0, t1, t2) = (1.0, t, t * t);

        let (u0, u1, u2) = (1.0, u, u * u);

        let mut dd = 0.0; // Result

        let v = zxyminmin;
        dd += v * t0 * u0;
        let v = -3.0 * zxminmin + 3.0 * zxminmax - 2.0 * zxyminmin - zxyminmax;
        dd += 2.0 * v * t0 * u1;
        let v = 2.0 * zxminmin - 2.0 * zxminmax + zxyminmin + zxyminmax;
        dd += 3.0 * v * t0 * u2;
        let v = -3.0 * zyminmin + 3.0 * zymaxmin - 2.0 * zxyminmin - zxymaxmin;
        dd += 2.0 * v * t1 * u0;
        #[rustfmt::skip]
        let v = 9.0 * zminmin - 9.0 * zmaxmin + 9.0 * zmaxmax - 9.0 * zminmax
            + 6.0 * zxminmin + 3.0 * zxmaxmin - 3.0 * zxmaxmax - 6.0 * zxminmax
            + 6.0 * zyminmin - 6.0 * zymaxmin - 3.0 * zymaxmax + 3.0 * zyminmax
            + 4.0 * zxyminmin + 2.0 * zxymaxmin + zxymaxmax + 2.0 * zxyminmax;
        dd += 4.0 * v * t1 * u1;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 4.0 * zxminmin - 2.0 * zxmaxmin + 2.0 * zxmaxmax + 4.0 * zxminmax
            - 3.0 * zyminmin + 3.0 * zymaxmin + 3.0 * zymaxmax - 3.0 * zyminmax
            - 2.0 * zxyminmin - zxymaxmin - zxymaxmax - 2.0 * zxyminmax;
        dd += 6.0 * v * t1 * u2;
        let v = 2.0 * zyminmin - 2.0 * zymaxmin + zxyminmin + zxymaxmin;
        dd += 3.0 * v * t2 * u0;
        #[rustfmt::skip]
        let v = -6.0 * zminmin + 6.0 * zmaxmin - 6.0 * zmaxmax + 6.0 * zminmax
            - 3.0 * zxminmin - 3.0 * zxmaxmin + 3.0 * zxmaxmax + 3.0 * zxminmax
            - 4.0 * zyminmin + 4.0 * zymaxmin + 2.0 * zymaxmax - 2.0 * zyminmax
            - 2.0 * zxyminmin - 2.0 * zxymaxmin - zxymaxmax - zxyminmax;
        dd += 6.0 * v * t2 * u1;
        #[rustfmt::skip]
        let v = 4.0 * zminmin - 4.0 * zmaxmin + 4.0 * zmaxmax - 4.0 * zminmax
            + 2.0 * zxminmin + 2.0 * zxmaxmin - 2.0 * zxmaxmax - 2.0 * zxminmax
            + 2.0 * zyminmin - 2.0 * zymaxmin - 2.0 * zymaxmax + 2.0 * zyminmax
            + zxyminmin + zxymaxmin + zxymaxmax + zxyminmax;
        dd += 9.0 * v * t2 * u2;
        dd = dd * dt * du;

        Ok(dd)
    }
}

/// Common calculation.
fn tu_cubic_values(x: f64, y: f64, xlo: f64, ylo: f64, dx: f64, dy: f64) -> (f64, f64, f64, f64) {
    let t = (x - xlo) / dx;
    let u = (y - ylo) / dy;
    let dt = dx.recip();
    let du = dy.recip();
    (t, u, dt, du)
}
