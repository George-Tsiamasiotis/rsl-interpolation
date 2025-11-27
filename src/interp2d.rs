use crate::types::check_if_inbounds;
use crate::{Accelerator, Cache, DynInterp2dType};
use crate::{DomainError, InterpolationError};

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
    type Interpolation2d: Interpolation2d<T> + Send + Sync;

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
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`, and the [`Cache`] `cache`.
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
    /// let mut cache = Cache::new();
    ///
    /// let z = interp.eval(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc, &mut cache)?;
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
        cache: &mut Cache<T>,
    ) -> Result<T, DomainError>
    where
        T: PartialOrd + Clone,
    {
        // Calculation is the same, with the added bounds check
        check_if_inbounds(xa, x.clone())?;
        check_if_inbounds(ya, y.clone())?;

        self.eval_extrap(xa, ya, za, x, y, xacc, yacc, cache)
    }

    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`, and the [`Cache`] `cache`.
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
    /// let mut cache = Cache::new();
    ///
    /// let z = interp.eval_extrap(&xa, &ya, &za, 3.0, 6.0, &mut xacc, &mut yacc, &mut cache)?;
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
        cache: &mut Cache<T>,
    ) -> Result<T, DomainError>;

    /// Returns the interpolated value `d = ∂z/∂x` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`, and the [`Cache`] `cache`.
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
    /// let mut cache = Cache::new();
    ///
    /// let dzdx = interp.eval_deriv_x(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc, &mut cache)?;
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
        cache: &mut Cache<T>,
    ) -> Result<T, DomainError>;

    /// Returns the interpolated value `d = ∂z/∂y` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`, and the [`Cache`] `cache`.
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
    /// let mut cache = Cache::new();
    ///
    /// let dzdy = interp.eval_deriv_y(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc, &mut cache)?;
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
        cache: &mut Cache<T>,
    ) -> Result<T, DomainError>;

    /// Returns the interpolated value `d = 𝜕²z/𝜕x²` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`, and the [`Cache`] `cache`.
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
    /// let mut cache = Cache::new();
    ///
    /// let dzdx2 = interp.eval_deriv_xx(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc, &mut cache)?;
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
        cache: &mut Cache<T>,
    ) -> Result<T, DomainError>;

    /// Returns the interpolated value `d = 𝜕²z/𝜕y²` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`, and the [`Cache`] `cache`.
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
    /// let mut cache = Cache::new();
    ///
    /// let dzdy2 = interp.eval_deriv_yy(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc, &mut cache)?;
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
        cache: &mut Cache<T>,
    ) -> Result<T, DomainError>;

    /// Returns the interpolated value `d = 𝜕²z/𝜕x𝜕y` for a given point (`x`, `y`), using the data arrays
    /// `xa`, `ya`, `za` and the [`Accelerators`] `xacc` and `yacc`, and the [`Cache`] `cache`.
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
    /// let mut cache = Cache::new();
    ///
    /// let dzdxy = interp.eval_deriv_xy(&xa, &ya, &za, 1.5, 3.0, &mut xacc, &mut yacc, &mut cache)?;
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
        cache: &mut Cache<T>,
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

/// Creates a [`DynInterp2dType`] of `typ` type.
///
/// Useful when `typ` is not known at compile time.
///
/// # Example
/// ```
/// # use rsl_interpolation::*;
/// #
/// # fn main() -> Result<(), InterpolationError> {
/// let xa = [0.0, 1.0, 2.0, 3.0];
/// let ya = [0.0, 2.0, 4.0, 6.0];
/// // z = x + y
/// let za = [
///     0.0, 1.0, 2.0, 3.0,
///     2.0, 3.0, 4.0, 5.0,
///     4.0, 5.0, 6.0, 7.0,
///     6.0, 7.0, 8.0, 9.0,
/// ];
/// let typ = "bicubic";
///
/// let interp2d_type = make_interp2d_type(typ)?;
/// let interp = interp2d_type.build(&xa, &ya, &za)?;
/// # Ok(())
/// # }
/// ```
pub fn make_interp2d_type<T>(typ: &str) -> Result<DynInterp2dType<T>, InterpolationError>
where
    T: crate::Num + ndarray_linalg::Lapack,
{
    use crate::*;

    match typ.to_lowercase().as_str() {
        "bilinear" => Ok(DynInterp2dType::new(Bilinear)),
        "bicubic" => Ok(DynInterp2dType::new(Bicubic)),
        _ => Err(InterpolationError::InvalidType(typ.into())),
    }
}

// ===============================================================================================

#[cfg(test)]
mod test {
    use super::*;
    use crate::*;

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

    #[test]
    fn test_z_get() {
        #[rustfmt::skip]
        let xa = [0.0, 1.0, 2.0];
        let ya = [0.0, 1.0, 2.0, 3.0];
        #[rustfmt::skip]
        let za = [
            0.0, 1.0, 2.0,
            3.0, 4.0, 5.0,
            6.0, 5.0, 4.0,
            3.0, 99.0, 1.0, // we want 99.0
        ];

        let (i, j) = (1, 3);
        let idx = z_get(&za, i, j, xa.len(), ya.len()).unwrap();
        let expected = 99.0;
        assert_eq!(idx, expected);
    }

    #[test]
    fn test_dyn_interp_type() {
        let xa = [0.0, 1.0, 2.0, 3.0];
        let ya = [0.0, 2.0, 4.0, 6.0];
        #[rustfmt::skip]
        let za = [
            0.0, 1.0, 2.0, 3.0,
            2.0, 3.0, 4.0, 5.0,
            4.0, 5.0, 6.0, 7.0,
            6.0, 7.0, 8.0, 9.0,
        ];
        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();
        let mut cache = Cache::new();

        let x = 0.5;
        let y = 1.0;
        let interp2d_type = DynInterp2dType::new(Bicubic);
        let interp2d = interp2d_type.build(&xa, &ya, &za).unwrap();
        interp2d
            .eval(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
    }
}
