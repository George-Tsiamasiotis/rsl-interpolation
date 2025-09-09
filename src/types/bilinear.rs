use std::marker::PhantomData;

use crate::Accelerator;
use crate::DomainError;
use crate::Interp2dType;
use crate::Interpolation2d;
use crate::InterpolationError;
use crate::interp2d::{acc_indices, partials, xy_grid_indices, z_grid_indices};
use crate::types::utils::check_if_inbounds;
use crate::types::utils::check2d_data;

const MIN_SIZE: usize = 2;

/// Bilinear Interpolation type.
///
/// The simplest type of 2d Interpolation.
#[doc(alias = "gsl_interp2d_bilinear")]
pub struct Bilinear;

impl<T> Interp2dType<T> for Bilinear
where
    T: crate::Num,
{
    type Interpolation2d = BilinearInterp<T>;

    /// Constructs a Bilinear Interpolator.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x + y, in column-major order
    /// let za = [
    ///     0.0, 1.0, 2.0,
    ///     2.0, 3.0, 4.0,
    ///     4.0, 5.0, 6.0,
    /// ];
    ///
    /// let interp = Bilinear.build(&xa, &ya, &za)?;
    /// # Ok(())
    /// # }
    /// ```
    fn build(&self, xa: &[T], ya: &[T], za: &[T]) -> Result<BilinearInterp<T>, InterpolationError> {
        check2d_data(xa, ya, za, MIN_SIZE)?;

        Ok(BilinearInterp {
            _variable_type: PhantomData,
        })
    }

    fn name(&self) -> &str {
        "Bilinear"
    }

    fn min_size(&self) -> usize {
        MIN_SIZE
    }
}

// ===============================================================================================

/// Bilinear Interpolator.
///
/// Provides all the evaluation methods.
///
/// Should be constructed through the [`Bilinear`] type.
pub struct BilinearInterp<T> {
    _variable_type: PhantomData<T>,
}

impl<T> Interpolation2d<T> for BilinearInterp<T>
where
    T: crate::Num,
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
        let (xi, yi) = acc_indices(xa, ya, x, y, xacc, yacc);
        let (xlo, xhi, ylo, yhi) = xy_grid_indices(xa, ya, xi, yi);
        let (zlolo, zlohi, zhilo, zhihi) = z_grid_indices(za, xa.len(), ya.len(), xi, yi)?;
        let (dx, dy) = partials(xlo, xhi, ylo, yhi);

        debug_assert!(dx > T::zero());
        debug_assert!(dy > T::zero());

        let t = (x - xlo) / dx;
        let u = (y - ylo) / dy;

        let one = T::one();
        let result = (one - t) * (one - u) * zlolo
            + t * (one - u) * zhilo
            + (one - t) * u * zlohi
            + t * u * zhihi;
        Ok(result)
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

        let (xi, yi) = acc_indices(xa, ya, x, y, xacc, yacc);
        let (xlo, xhi, ylo, yhi) = xy_grid_indices(xa, ya, xi, yi);
        let (zlolo, zlohi, zhilo, zhihi) = z_grid_indices(za, xa.len(), ya.len(), xi, yi)?;
        let (dx, dy) = partials(xlo, xhi, ylo, yhi);

        debug_assert!(dx > T::zero());
        debug_assert!(dy > T::zero());

        let one = T::one();
        let dt = one / dx;
        let u = (y - ylo) / dy;

        let result = dt * (-(one - u) * zlolo + (one - u) * zhilo - u * zlohi + u * zhihi);
        Ok(result)
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

        let (xi, yi) = acc_indices(xa, ya, x, y, xacc, yacc);
        let (xlo, xhi, ylo, yhi) = xy_grid_indices(xa, ya, xi, yi);
        let (zlolo, zlohi, zhilo, zhihi) = z_grid_indices(za, xa.len(), ya.len(), xi, yi)?;
        let (dx, dy) = partials(xlo, xhi, ylo, yhi);

        debug_assert!(dx > T::zero());
        debug_assert!(dy > T::zero());

        let one = T::one();
        let t = (x - xlo) / dx;
        let du = one / dy;

        let result = du * (-(one - t) * zlolo - t * zhilo + (one - t) * zlohi + t * zhihi);
        Ok(result)
    }

    #[allow(unused_variables)]
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

        Ok(T::zero())
    }

    #[allow(unused_variables)]
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

        Ok(T::zero())
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

        let (xi, yi) = acc_indices(xa, ya, x, y, xacc, yacc);
        let (xlo, xhi, ylo, yhi) = xy_grid_indices(xa, ya, xi, yi);
        let (zlolo, zlohi, zhilo, zhihi) = z_grid_indices(za, xa.len(), ya.len(), xi, yi)?;
        let (dx, dy) = partials(xlo, xhi, ylo, yhi);

        let one = T::one();
        let dt = one / dx;
        let du = one / dy;

        let result = dt * du * (zlolo - zhilo - zlohi + zhihi);
        Ok(result)
    }
}
