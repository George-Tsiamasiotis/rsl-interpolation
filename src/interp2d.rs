use crate::types::check_if_inbounds;
use crate::{Accelerator, DomainError, InterpolationError};

/// Defines the required methods for every 2d Interpolation type.
pub trait Interpolation2d<T>
where
    T: num::Float + std::fmt::Debug,
{
    /// The minimum number of points required by the interpolator. For example, bicubic
    /// interpolation requires a minimum of 4 points.
    const MIN_SIZE: usize;

    /// The name of the interpolator.
    const NAME: &'static str;

    /// Creates a new 2d Interpolator for the data (xa, ya, za), where xa and ya are slices of the
    /// x and y grid points and za is an array of functions values of len(xa)*len(ya).
    #[doc(alias = "gsl_interp2d_init")]
    fn new(xa: &[T], ya: &[T], za: &[T]) -> Result<Self, InterpolationError>
    where
        Self: Sized;

    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the Accelerators `xacc` and `yacc`.
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`DomainError`]: struct.DomainError.html
    #[doc(alias = "gsl_interp2d_eval")]
    #[doc(alias = "gsl_interp2d_eval_e")]
    #[allow(clippy::too_many_arguments)]
    fn eval(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        // Calculation is the same, with the added bounds check
        check_if_inbounds(xa, x)?;
        check_if_inbounds(ya, y)?;

        self.eval_extrap(xa, ya, za, x, y, xacc, yacc)
    }

    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the Accelerators `xacc` and `yacc`.
    ///
    /// # Note
    ///
    /// This function performs *no bound checking*,  se when `x` is outside the range of `xa` or y
    /// is outside the range of `ya`, extrapolation is performed.
    #[doc(alias = "gsl_interp2d_eval_extrap")]
    #[doc(alias = "gsl_interp2d_eval_extrap_e")]
    #[allow(clippy::too_many_arguments)]
    fn eval_extrap(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError>;

    /// Returns the interpolated value `d = ∂z/∂x` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the Accelerators `xacc` and `yacc`.
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`DomainError`]: struct.DomainError.html
    #[doc(alias = "gsl_interp2d_eval_deriv_x")]
    #[doc(alias = "gsl_interp2d_eval_deriv_x_e")]
    #[allow(clippy::too_many_arguments)]
    fn eval_deriv_x(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError>;

    /// Returns the interpolated value `d = ∂z/∂y` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the Accelerators `xacc` and `yacc`.
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`DomainError`]: struct.DomainError.html
    #[doc(alias = "gsl_interp2d_eval_deriv_y")]
    #[doc(alias = "gsl_interp2d_eval_deriv_y_e")]
    #[allow(clippy::too_many_arguments)]
    fn eval_deriv_y(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError>;

    /// Returns the interpolated value `d = 𝜕²z/𝜕x²` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the Accelerators `xacc` and `yacc`.
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`DomainError`]: struct.DomainError.html
    #[doc(alias = "gsl_interp2d_eval_deriv_xx")]
    #[doc(alias = "gsl_interp2d_eval_deriv_xx_e")]
    #[allow(clippy::too_many_arguments)]
    fn eval_deriv_xx(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError>;

    /// Returns the interpolated value `d = 𝜕²z/𝜕y²` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the Accelerators `xacc` and `yacc`.
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`DomainError`]: struct.DomainError.html
    #[doc(alias = "gsl_interp2d_eval_deriv_yy")]
    #[doc(alias = "gsl_interp2d_eval_deriv_yy_e")]
    #[allow(clippy::too_many_arguments)]
    fn eval_deriv_yy(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError>;

    /// Returns the interpolated value `d = 𝜕²z/𝜕x𝜕y` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the Accelerators `xacc` and `yacc`.
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`DomainError`]: struct.DomainError.html
    #[doc(alias = "gsl_interp2d_eval_deriv_xy")]
    #[doc(alias = "gsl_interp2d_eval_deriv_xy_e")]
    #[allow(clippy::too_many_arguments)]
    fn eval_deriv_xy(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError>;

    // TODO: decide whether to make these trait methods or not

    /// Sets the value `z` of grid point (`i`, `j`) of the array `za` to `z`.
    #[doc(alias = "gsl_inter2d_set")]
    #[allow(unused_variables)] // FIXME:
    fn set(za: &mut [T], z: T, i: usize, j: usize) -> Result<(), InterpolationError> {
        todo!()
    }

    /// Returns the value `z` of grid point (`i`, `j`) of the array `za` to `z`.
    #[doc(alias = "gsl_inter2d_get")]
    #[allow(unused_variables)] // FIXME:
    fn get(za: &mut [T], z: T, i: usize, j: usize) -> Result<(), InterpolationError> {
        todo!()
    }

    /// Returns the index corresponding to the grid point (`i`, `j`). The index is given by
    /// `j*len(x) + i`
    #[doc(alias = "gsl_interp2d_idx")]
    fn idx(i: usize, j: usize, ilen: usize, jlen: usize) -> Result<usize, DomainError> {
        idx(i, j, ilen, jlen)
    }
}

pub(crate) fn idx(xi: usize, yi: usize, xlen: usize, ylen: usize) -> Result<usize, DomainError> {
    if (xi >= xlen) | (yi >= ylen) {
        Err(DomainError)
    } else {
        Ok(yi * xlen + xi)
    }
}

/// Common calculation to evaluation functions
pub(crate) fn acc_indeces<T>(
    xa: &[T],
    ya: &[T],
    x: T,
    y: T,
    xacc: &mut Accelerator,
    yacc: &mut Accelerator,
) -> (usize, usize)
where
    T: num::Float + std::fmt::Debug,
{
    let xi = xacc.find(xa, x);
    let yi = yacc.find(ya, y);
    (xi, yi)
}

/// Common calculation to evaluation functions
pub(crate) fn xy_grid_indeces<T>(xa: &[T], ya: &[T], xi: usize, yi: usize) -> (T, T, T, T)
where
    T: num::Float + std::fmt::Debug,
{
    let xlo = xa[xi];
    let xhi = xa[xi + 1];
    let ylo = ya[yi];
    let yhi = ya[yi + 1];
    (xlo, xhi, ylo, yhi)
}

/// Common calculation to evaluation functions
pub(crate) fn z_grid_indeces<T>(
    za: &[T],
    xlen: usize,
    ylen: usize,
    xi: usize,
    yi: usize,
) -> Result<(T, T, T, T), DomainError>
where
    T: num::Float + std::fmt::Debug,
{
    let zlolo = za[idx(xi, yi, xlen, ylen)?];
    let zlohi = za[idx(xi, yi + 1, xlen, ylen)?];
    let zhilo = za[idx(xi + 1, yi, xlen, ylen)?];
    let zhihi = za[idx(xi + 1, yi + 1, xlen, ylen)?];
    Ok((zlolo, zlohi, zhilo, zhihi))
}

/// Common calculation to evaluation functions
pub(crate) fn partials<T>(xlo: T, xhi: T, ylo: T, yhi: T) -> (T, T)
where
    T: num::Float + std::fmt::Debug,
{
    let dx = xhi - xlo;
    let dy = yhi - ylo;
    (dx, dy)
}
