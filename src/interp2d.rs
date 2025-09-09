use crate::types::check_if_inbounds;
use crate::{Accelerator, DomainError, InterpolationError};

/// Representation of a 2D Interpolation Type.
///
/// > # **Important**
/// >
/// > The `za` array must be defined in **column-major (Fortran)** style. This is done to comply
/// > with GSL's interface.
/// >
///
/// For 2d interpolation, 2 separate [`Accelerators`] are required for each of the grid variables.
///
/// [`Accelerators`]: Accelerator
pub trait Interp2dType<T> {
    /// The returned 2D Interpolator, containing the calculated coefficients and providing the
    /// evaluation methods.
    type Interpolation2d: Interpolation2d<T>;

    /// Creates a 2D Interpolator from the data arrays `xa`, `ya` and `za`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0, 3.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0];
    /// // z = x + y
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
    fn build(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
    ) -> Result<Self::Interpolation2d, InterpolationError>;

    /// Returns the name of the Interpolator.
    #[doc(alias = "gsl_interp_name")]
    fn name(&self) -> &str;

    /// Returns the minimum number of points required by the Interpolator.
    #[doc(alias = "gsl_interp_min_size")]
    fn min_size(&self) -> usize;
}

/// Defines the required evaluation methods.
#[allow(clippy::too_many_arguments)]
pub trait Interpolation2d<T> {
    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Note
    ///
    /// This function only performs the bounds check, and then calls `eval_extrap()`, where the
    /// actual evaluation is implemented.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x + y
    /// let za = [
    ///     0.0, 1.0, 2.0,
    ///     2.0, 3.0, 4.0,
    ///     4.0, 5.0, 6.0,
    /// ];
    /// let interp = Bilinear.build(&xa, &ya, &za)?;
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let z = interp.eval(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc)?;
    ///
    /// assert_eq!(z, 4.5);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`Accelerators`]: Accelerator
    #[doc(alias = "gsl_interp2d_eval")]
    #[doc(alias = "gsl_interp2d_eval_e")]
    fn eval(
        &self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError>
    where
        T: PartialOrd + Clone,
    {
        // Calculation is the same, with the added bounds check
        check_if_inbounds(xa, x.clone())?;
        check_if_inbounds(ya, y.clone())?;

        self.eval_extrap(xa, ya, za, x, y, xacc, yacc)
    }

    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerators`]` `xacc` and `yacc`.
    ///
    /// # Note
    ///
    /// This function performs *no bound checking*, so when `x` is outside the range of `xa` or y
    /// is outside the range of `ya`, extrapolation is performed.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x + y
    /// let za = [
    ///     0.0, 1.0, 2.0,
    ///     2.0, 3.0, 4.0,
    ///     4.0, 5.0, 6.0,
    /// ];
    /// let interp = Bilinear.build(&xa, &ya, &za)?;
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let z = interp.eval_extrap(&xa, &ya, &za, 3.0, 6.0, &mut xacc, &mut yacc)?;
    ///
    /// assert_eq!(z, 9.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// [`Accelerators`]: Accelerator
    #[doc(alias = "gsl_interp2d_eval_extrap")]
    #[doc(alias = "gsl_interp2d_eval_extrap_e")]
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
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let interp = Bilinear.build(&xa, &ya, &za)?;
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let dzdx = interp.eval_deriv_x(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc)?;
    ///
    /// assert_eq!(dzdx, 3.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`Accelerators`]: Accelerator
    #[doc(alias = "gsl_interp2d_eval_deriv_x")]
    #[doc(alias = "gsl_interp2d_eval_deriv_x_e")]
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
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let interp = Bilinear.build(&xa, &ya, &za)?;
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let dzdy = interp.eval_deriv_y(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc)?;
    ///
    /// assert_eq!(dzdy, 6.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`Accelerators`]: Accelerator
    #[doc(alias = "gsl_interp2d_eval_deriv_y")]
    #[doc(alias = "gsl_interp2d_eval_deriv_y_e")]
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
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let interp = Bilinear.build(&xa, &ya, &za)?;
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let dzdx2 = interp.eval_deriv_xx(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc)?;
    ///
    /// assert_eq!(dzdx2, 0.0); // Linear Interpolation!
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`Accelerators`]: Accelerator
    #[doc(alias = "gsl_interp2d_eval_deriv_xx")]
    #[doc(alias = "gsl_interp2d_eval_deriv_xx_e")]
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
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let interp = Bilinear.build(&xa, &ya, &za)?;
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let dzdy2 = interp.eval_deriv_yy(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc)?;
    ///
    /// assert_eq!(dzdy2, 0.0); // Linear Interpolation!
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`Accelerators`]: Accelerator
    #[doc(alias = "gsl_interp2d_eval_deriv_yy")]
    #[doc(alias = "gsl_interp2d_eval_deriv_yy_e")]
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
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let interp = Bilinear.build(&xa, &ya, &za)?;
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let dzdxy = interp.eval_deriv_xy(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc)?;
    ///
    /// assert_eq!(dzdxy, 0.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`DomainError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    ///
    /// [`Accelerators`]: Accelerator
    #[doc(alias = "gsl_interp2d_eval_deriv_xy")]
    #[doc(alias = "gsl_interp2d_eval_deriv_xy_e")]
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
}

/// Returns the index corresponding to the grid point (`i`, `j`). The index is given by
/// `j*len(x) + i`
///
/// > # Important
/// >
/// > The `za` array is indexed in column-major style (Fortran style), so it must be defined
/// > accordingly.
///
/// # Example
///
/// ```
/// # use rsl_interpolation::*;
/// #
/// # fn main() -> Result<(), DomainError>{
/// let xa = [0.0, 1.0];
/// let ya = [0.0, 2.0];
/// let za = [
///     0.0, 1.0, // <- This one
///     2.0, 3.0,
/// ];
/// let za_index = z_idx(0, 1, xa.len(), ya.len())?;
/// assert_eq!(za_index, 2);
/// #    Ok(())
/// }
#[doc(alias = "gsl_interp2d_idx")]
pub fn z_idx(xi: usize, yi: usize, xlen: usize, ylen: usize) -> Result<usize, DomainError> {
    if (xi >= xlen) | (yi >= ylen) {
        Err(DomainError)
    } else {
        Ok(yi * xlen + xi)
    }
}

/// Sets the value `z` of grid point (`i`, `j`) of the array `za` to `z`.
///
/// > # Important
/// >
/// > The `za` array is indexed in column-major style (Fortran style), so it must be defined
/// > accordingly.
///
/// # Example
///
/// ```
/// # use rsl_interpolation::*;
/// #
/// # fn main() -> Result<(), DomainError>{
/// let xa = [0.0, 1.0];
/// let ya = [0.0, 2.0];
/// let mut za = [
///     0.0, 1.0, // <- We set this one
///     2.0, 3.0,
/// ];
/// z_set(&mut za, 10.0, 0, 1, xa.len(), ya.len())?;
/// assert_eq!(za[2], 10.0);
/// #    Ok(())
/// }
#[doc(alias = "gsl_inter2d_set")]
pub fn z_set<T>(
    za: &mut [T],
    z: T,
    i: usize,
    j: usize,
    xlen: usize,
    ylen: usize,
) -> Result<(), DomainError>
where
    T: crate::Num,
{
    if (i >= xlen) | (j >= ylen) {
        return Err(DomainError);
    };

    za[z_idx(i, j, xlen, ylen)?] = z;

    Ok(())
}

/// Returns the value `z` of grid point (`i`, `j`) of the array `za` to `z`.
///
/// > # Important
/// >
/// > The `za` array is indexed in column-major style (Fortran style), so it must be defined
/// > accordingly.
///
/// # Example
///
/// ```
/// # use rsl_interpolation::*;
/// #
/// # fn main() -> Result<(), DomainError>{
/// let xa = [0.0, 1.0];
/// let ya = [0.0, 2.0];
/// let za = [
///     0.0, 10.0, // <- We want this one
///     2.0, 3.0,
/// ];
/// let g = z_get(&za, 1, 0, xa.len(), ya.len())?;
/// assert_eq!(g, 10.0);
/// #    Ok(())
/// }
#[doc(alias = "gsl_inter2d_get")]
pub fn z_get<T>(za: &[T], i: usize, j: usize, xlen: usize, ylen: usize) -> Result<T, DomainError>
where
    T: crate::Num,
{
    if (i >= xlen) | (j >= ylen) {
        return Err(DomainError);
    };

    Ok(za[z_idx(i, j, xlen, ylen)?])
}

// ===============================================================================================

/// Common calculation to evaluation functions
pub(crate) fn acc_indices<T>(
    xa: &[T],
    ya: &[T],
    x: T,
    y: T,
    xacc: &mut Accelerator,
    yacc: &mut Accelerator,
) -> (usize, usize)
where
    T: crate::Num,
{
    let xi = xacc.find(xa, x);
    let yi = yacc.find(ya, y);
    (xi, yi)
}

/// Common calculation to evaluation functions
pub(crate) fn xy_grid_indices<T>(xa: &[T], ya: &[T], xi: usize, yi: usize) -> (T, T, T, T)
where
    T: crate::Num,
{
    let xlo = xa[xi];
    let xhi = xa[xi + 1];
    let ylo = ya[yi];
    let yhi = ya[yi + 1];
    (xlo, xhi, ylo, yhi)
}

/// Common calculation to evaluation functions
pub(crate) fn z_grid_indices<T>(
    za: &[T],
    xlen: usize,
    ylen: usize,
    xi: usize,
    yi: usize,
) -> Result<(T, T, T, T), DomainError>
where
    T: crate::Num,
{
    let zlolo = za[z_idx(xi, yi, xlen, ylen)?];
    let zlohi = za[z_idx(xi, yi + 1, xlen, ylen)?];
    let zhilo = za[z_idx(xi + 1, yi, xlen, ylen)?];
    let zhihi = za[z_idx(xi + 1, yi + 1, xlen, ylen)?];
    Ok((zlolo, zlohi, zhilo, zhihi))
}

/// Common calculation to evaluation functions
pub(crate) fn partials<T>(xlo: T, xhi: T, ylo: T, yhi: T) -> (T, T)
where
    T: crate::Num,
{
    let dx = xhi - xlo;
    let dy = yhi - ylo;
    (dx, dy)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_set() {
        let xa = [0.0, 1.0];
        let ya = [0.0, 2.0];

        #[rustfmt::skip]
        let mut za = [
            0.0, 1.0,
            1.0, 0.5,
        ];

        let za00 = 100.0;
        let za01 = 300.0;
        let za10 = 200.0;
        let za11 = 400.0;

        let xlen = xa.len();
        let ylen = ya.len();

        z_set(&mut za, za00, 0, 0, xlen, ylen).unwrap();
        z_set(&mut za, za01, 0, 1, xlen, ylen).unwrap();
        z_set(&mut za, za10, 1, 0, xlen, ylen).unwrap();
        z_set(&mut za, za11, 1, 1, xlen, ylen).unwrap();

        assert_eq!(za, [100.0, 200.0, 300.0, 400.0,])
    }
}
