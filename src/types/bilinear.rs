//! Definition of Bilinear Interpolator.

use crate::Accelerator2d;
use crate::{Domain2dError, InterpolatorError};
use crate::{Interpolation2d, Interpolator2d};
use crate::{check_if_inbounds2d, check2d_data};

const MIN_SIZE: usize = 2;

/// Bilinear Interpolator.
///
/// This interpolation method does not require any additional memory.
#[doc(alias = "gsl_interp2d_bilinear")]
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct BilinearInterpolator;

impl Interpolator2d for BilinearInterpolator {
    /// Constructs a Bilinear Interpolator.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x + y, in column-major order
    /// let za = [
    ///     0.0, 1.0, 2.0,
    ///     2.0, 3.0, 4.0,
    ///     4.0, 5.0, 6.0,
    /// ];
    ///
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// - [`InterpolatorError::UnsortedDataset`]: `xa` or `ya` are not monotonically increasing.
    /// - [`InterpolatorError::NotEnoughPoints`]: length of `xa` or `ya` is less that 2.
    /// - [`InterpolatorError::ZGridMismatch`]: `xa.len()*ya.len() != za.len()`.
    fn build(xa: &[f64], ya: &[f64], za: &[f64]) -> Result<Self, InterpolatorError> {
        check2d_data(xa, ya, za, MIN_SIZE)?;
        Ok(Self)
    }

    fn name(&self) -> &str {
        "Bilinear"
    }

    fn min_size(&self) -> usize {
        MIN_SIZE
    }
}

// ===============================================================================================

impl Interpolation2d for BilinearInterpolator {
    fn eval_extrap(
        &self,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        let is_uptodate = acc.is_uptodate(xa, ya, x, y);
        if !is_uptodate {
            acc.update_step1(xa, ya, za);
        }

        let (xlo, _, ylo, _) = acc.get_xy_grid_values();
        let (zminmin, zminmax, zmaxmin, zmaxmax) = acc.get_z_grid_values();
        let (dx, dy) = acc.get_partials();

        debug_assert!(dx > 0.0);
        debug_assert!(dy > 0.0);

        let t = (x - xlo) / dx;
        let u = (y - ylo) / dy;

        let result = (1.0 - t) * (1.0 - u) * zminmin
            + t * (1.0 - u) * zmaxmin
            + (1.0 - t) * u * zminmax
            + t * u * zmaxmax;
        Ok(result)
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

        let (_, _, ylo, _) = acc.get_xy_grid_values();
        let (zminmin, zminmax, zmaxmin, zmaxmax) = acc.get_z_grid_values();
        let (dx, dy) = acc.get_partials();

        debug_assert!(dx > 0.0);
        debug_assert!(dy > 0.0);

        let dt = 1.0 / dx;
        let u = (y - ylo) / dy;

        let result = dt * (-(1.0 - u) * zminmin + (1.0 - u) * zmaxmin - u * zminmax + u * zmaxmax);
        Ok(result)
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

        let (xlo, _, _, _) = acc.get_xy_grid_values();
        let (zminmin, zminmax, zmaxmin, zmaxmax) = acc.get_z_grid_values();
        let (dx, dy) = acc.get_partials();

        debug_assert!(dx > 0.0);
        debug_assert!(dy > 0.0);

        let t = (x - xlo) / dx;
        let du = 1.0 / dy;

        let result = du * (-(1.0 - t) * zminmin - t * zmaxmin + (1.0 - t) * zminmax + t * zmaxmax);
        Ok(result)
    }

    fn eval_deriv_xx(
        &self,
        xa: &[f64],
        ya: &[f64],
        _: &[f64],
        x: f64,
        y: f64,
        _: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        check_if_inbounds2d(xa, ya, x, y)?;
        Ok(0.0)
    }

    fn eval_deriv_yy(
        &self,
        xa: &[f64],
        ya: &[f64],
        _: &[f64],
        x: f64,
        y: f64,
        _: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        check_if_inbounds2d(xa, ya, x, y)?;

        Ok(0.0)
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

        let (zminmin, zminmax, zmaxmin, zmaxmax) = acc.get_z_grid_values();
        let (dx, dy) = acc.get_partials();

        let dt = 1.0 / dx;
        let du = 1.0 / dy;

        let result = dt * du * (zminmin - zmaxmin - zminmax + zmaxmax);
        Ok(result)
    }
}
