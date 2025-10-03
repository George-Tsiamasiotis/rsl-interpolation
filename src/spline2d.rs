use crate::Accelerator;
use crate::DynInterp2dType;
use crate::Interp2dType;
use crate::Interpolation2d;
use crate::{DomainError, InterpolationError};

/// 2D Higher level interface.
///
/// A 2D Spline owns the data it is constructed with, and provides the same evaluation methods as the
/// lower-level Interpolator object, without needing to provide the data arrays in every call.
///
/// # Example
/// ```
/// # use rsl_interpolation::*;
/// #
/// # fn main() -> Result<(), InterpolationError>{
/// let mut xacc = Accelerator::new();
/// let mut yacc = Accelerator::new();
///
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
///
/// let typ = Bicubic;
/// let spline = Spline2d::new(typ, &xa, &ya, &za)?;
///
/// let (x, y) = (2.5, 4.1);
/// let y_interp = interp.eval(&xa, &ya, &za, x, y, &mut xacc, &mut yacc)?;
/// let y_spline = spline.eval(x, y, &mut xacc, &mut yacc)?;
///
/// assert_eq!(y_interp, y_spline);
/// #
/// # Ok(())
/// # }
/// ```
pub struct Spline2d<I, T>
where
    I: Interp2dType<T>,
{
    /// The lower-level [`2D Interpolator`].
    ///
    /// [`2D Interpolator`]: Interpolation2d#implementors
    pub interp: I::Interpolation2d,
    /// The owned x data.
    pub xa: Box<[T]>,
    /// The owned y data.
    pub ya: Box<[T]>,
    /// The owned z data.
    pub za: Box<[T]>,
    name: Box<str>,
    min_size: usize,
}

impl<I, T> Spline2d<I, T>
where
    I: Interp2dType<T>,
{
    /// Constructs a 2D Spline of a 2D Interpolation type `typ` from the data arrays `xa`, `ya` and
    /// `za`.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
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
    /// let typ = Bicubic;
    /// let spline = Spline2d::new(typ, &xa, &ya, &za)?;
    /// #
    /// # Ok(())
    /// # }
    /// ```
    #[doc(alias = "gsl_spline2d_init")]
    pub fn new(typ: I, xa: &[T], ya: &[T], za: &[T]) -> Result<Self, InterpolationError>
    where
        T: Clone,
    {
        let spline = Self {
            interp: typ.build(xa, ya, za)?,
            xa: xa.into(),
            ya: ya.into(),
            za: za.into(),
            name: typ.name().into(),
            min_size: typ.min_size(),
        };

        Ok(spline)
    }

    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the
    /// [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x + y
    /// let za = [
    ///     0.0, 1.0, 2.0,
    ///     2.0, 3.0, 4.0,
    ///     4.0, 5.0, 6.0,
    /// ];
    ///
    /// let typ = Bilinear;
    /// let spline = Spline2d::new(typ, &xa, &ya, &za)?;
    ///
    /// let z = spline.eval(1.5, 3.0, &mut xacc, &mut yacc)?;
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
    #[doc(alias = "gsl_spline2d_eval")]
    #[doc(alias = "gsl_spline2d_eval_e")]
    pub fn eval(
        &self,
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError>
    where
        T: crate::Num,
    {
        self.interp
            .eval(&self.xa, &self.ya, &self.za, x, y, xacc, yacc)
    }

    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the
    /// [`Accelerators`]` `xacc` and `yacc`.
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
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x + y
    /// let za = [
    ///     0.0, 1.0, 2.0,
    ///     2.0, 3.0, 4.0,
    ///     4.0, 5.0, 6.0,
    /// ];
    ///
    /// let typ = Bilinear;
    /// let spline = Spline2d::new(typ, &xa, &ya, &za)?;
    ///
    /// let z = spline.eval_extrap(3.0, 6.0, &mut xacc, &mut yacc)?;
    ///
    /// assert_eq!(z, 9.0);
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
    #[doc(alias = "gsl_spline2d_eval_extrap")]
    #[doc(alias = "gsl_spline2d_eval_extrap_e")]
    pub fn eval_extrap(
        &self,
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        self.interp
            .eval_extrap(&self.xa, &self.ya, &self.za, x, y, xacc, yacc)
    }

    /// Returns the interpolated value `d = ∂z/∂x` for a given point (`x`, `y`), using the
    /// [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    ///
    /// let typ = Bilinear;
    /// let spline = Spline2d::new(typ, &xa, &ya, &za)?;
    ///
    /// let dzdx = spline.eval_deriv_x(1.5, 3.0, &mut xacc, &mut yacc)?;
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
    #[doc(alias = "gsl_spline2d_eval_deriv_x")]
    #[doc(alias = "gsl_spline2d_eval_deriv_x_e")]
    pub fn eval_deriv_x(
        &self,
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        self.interp
            .eval_deriv_x(&self.xa, &self.ya, &self.za, x, y, xacc, yacc)
    }

    /// Returns the interpolated value `d = ∂z/∂y` for a given point (`x`, `y`), using the
    /// [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    ///
    /// let typ = Bilinear;
    /// let spline = Spline2d::new(typ, &xa, &ya, &za)?;
    ///
    /// let dzdx = spline.eval_deriv_y(1.5, 3.0, &mut xacc, &mut yacc)?;
    ///
    /// assert_eq!(dzdx, 6.0);
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
    #[doc(alias = "gsl_spline2d_eval_deriv_y")]
    #[doc(alias = "gsl_spline2d_eval_deriv_y_e")]
    pub fn eval_deriv_y(
        &self,
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        self.interp
            .eval_deriv_y(&self.xa, &self.ya, &self.za, x, y, xacc, yacc)
    }

    /// Returns the interpolated value `d = 𝜕²z/𝜕x²` for a given point (`x`, `y`), using the
    /// [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    ///
    /// let typ = Bilinear;
    /// let spline = Spline2d::new(typ, &xa, &ya, &za)?;
    ///
    /// let dzdx2 = spline.eval_deriv_xx(1.5, 3.0, &mut xacc, &mut yacc)?;
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
    pub fn eval_deriv_xx(
        &self,
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        self.interp
            .eval_deriv_xx(&self.xa, &self.ya, &self.za, x, y, xacc, yacc)
    }

    /// Returns the interpolated value `d = 𝜕²z/𝜕x²` for a given point (`x`, `y`), using the
    /// [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    /// let typ = Bilinear;
    /// let spline = Spline2d::new(typ, &xa, &ya, &za)?;
    ///
    /// let dzdx2 = spline.eval_deriv_yy(1.5, 3.0, &mut xacc, &mut yacc)?;
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
    #[doc(alias = "gsl_interp2d_eval_deriv_yy")]
    #[doc(alias = "gsl_interp2d_eval_deriv_yy_e")]
    pub fn eval_deriv_yy(
        &self,
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        self.interp
            .eval_deriv_yy(&self.xa, &self.ya, &self.za, x, y, xacc, yacc)
    }

    /// Returns the interpolated value `d = 𝜕²z/𝜕x𝜕y` for a given point (`x`, `y`), using the
    /// [`Accelerators`] `xacc` and `yacc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError>{
    /// let mut xacc = Accelerator::new();
    /// let mut yacc = Accelerator::new();
    ///
    /// let xa = [0.0, 1.0, 2.0];
    /// let ya = [0.0, 2.0, 4.0];
    /// // z = x² + y²
    /// let za = [
    ///      0.0,  1.0,  4.0,
    ///      4.0,  5.0,  8.0,
    ///     16.0, 17.0, 20.0,
    /// ];
    ///
    /// let typ = Bilinear;
    /// let spline = Spline2d::new(typ, &xa, &ya, &za)?;
    ///
    /// let dzdxy = spline.eval_deriv_xy(1.5, 3.0, &mut xacc, &mut yacc)?;
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
    pub fn eval_deriv_xy(
        &self,
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<T, DomainError> {
        self.interp
            .eval_deriv_xy(&self.xa, &self.ya, &self.za, x, y, xacc, yacc)
    }

    /// Returns the name of the Interpolator.
    #[doc(alias = "gsl_interp2d_name")]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the minimum number of points required by the Interpolator.
    #[doc(alias = "gsl_interp2d_min_size")]
    pub fn min_size(&self) -> usize {
        self.min_size
    }
}

/// 2D Spline with runtime-determined Interpolation Type.
pub type DynSpline2d<T> = Spline2d<DynInterp2dType<T>, T>;

impl<T> DynSpline2d<T> {
    /// Constructs a 2d Spline of a dynamic 2d Interpolation type `typ` from the data arrays `xa`,
    /// `ya` and `za`.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolationError> {
    /// let xa = [0.0, 1.0, 2.0, 3.0];
    /// let ya = [0.0, 2.0, 4.0, 6.0];
    /// // z = x + y, in column-major order
    /// let za = [
    ///     0.0, 1.0, 2.0, 3.0,
    ///     2.0, 3.0, 4.0, 5.0,
    ///     4.0, 5.0, 6.0, 7.0,
    ///     6.0, 7.0, 8.0, 9.0,
    /// ];
    /// let typ = "bicubic";
    ///
    /// let spline = match typ {
    ///     "bilinear" => Spline2d::new_dyn(Bilinear, &xa, &ya, &za)?,
    ///     "bicubic" => Spline2d::new_dyn(Bicubic, &xa, &ya, &za)?,
    ///     _ => unreachable!()
    /// };
    /// # Ok(())
    /// # }
    /// ```
    #[doc(alias = "gsl_spline2d_init")]
    pub fn new_dyn<I>(typ: I, xa: &[T], ya: &[T], za: &[T]) -> Result<Self, InterpolationError>
    where
        T: Clone,
        I: Interp2dType<T> + 'static,
        I::Interpolation2d: 'static,
    {
        Self::new(DynInterp2dType::new(typ), xa, ya, za)
    }
}

/// Creates a [`DynSpline2d`] of `typ` type.
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
/// // z = x + y, in column-major order
/// let za = [
///     0.0, 1.0, 2.0, 3.0,
///     2.0, 3.0, 4.0, 5.0,
///     4.0, 5.0, 6.0, 7.0,
///     6.0, 7.0, 8.0, 9.0,
/// ];
/// let typ = "bicubic";
///
/// let spline = make_spline2d(typ, &xa, &ya, &za)?;
/// # Ok(())
/// # }
/// ```
pub fn make_spline2d<T>(
    typ: &str,
    xa: &[T],
    ya: &[T],
    za: &[T],
) -> Result<DynSpline2d<T>, InterpolationError>
where
    T: crate::Num + ndarray_linalg::Lapack,
{
    use crate::*;

    match typ.to_lowercase().as_str() {
        "bilinear" => Ok(DynSpline2d::new_dyn(Bilinear, xa, ya, za)?),
        "bicubic" => Ok(DynSpline2d::new_dyn(Bicubic, xa, ya, za)?),
        _ => Err(InterpolationError::InvalidType(typ.into())),
    }
}

#[cfg(test)]
mod test {
    use crate::tests::build_comparator;
    use crate::*;

    #[test]
    fn test_spline2d_creation() {
        let xa = [0.0, 1.0, 2.0];
        let ya = [0.0, 2.0, 4.0];
        let za = [0.0, 1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 5.0, 6.0];

        let spline2d = Spline2d::new(Bilinear, &xa, &ya, &za).unwrap();
        let _: &str = spline2d.name();
        let _: usize = spline2d.min_size();
    }

    #[test]
    fn test_spline2d_eval() {
        let comp = build_comparator::<f64>();

        let mut xacc = Accelerator::new();
        let mut yacc = Accelerator::new();

        let xa = [0.0, 1.0, 2.0, 3.0];
        let ya = [0.0, 1.0, 2.0, 3.0];
        #[rustfmt::skip]
        let za = [
            1.0, 1.1, 1.2, 1.3,
            1.1, 1.2, 1.3, 1.4,
            1.2, 1.3, 1.4, 1.5,
            1.3, 1.4, 1.5, 1.6,
        ];

        let spline2d = Spline2d::new(Bicubic, &xa, &ya, &za).unwrap();

        let (x, y) = (1.5, 1.5);
        let z = spline2d.eval(x, y, &mut xacc, &mut yacc).unwrap();
        let dzdx = spline2d.eval_deriv_x(x, y, &mut xacc, &mut yacc).unwrap();
        let dzdy = spline2d.eval_deriv_y(x, y, &mut xacc, &mut yacc).unwrap();
        let dzdx2 = spline2d.eval_deriv_xx(x, y, &mut xacc, &mut yacc).unwrap();
        let dzdy2 = spline2d.eval_deriv_yy(x, y, &mut xacc, &mut yacc).unwrap();
        let dzdxy = spline2d.eval_deriv_xy(x, y, &mut xacc, &mut yacc).unwrap();

        assert!(comp.is_close(z, 1.3));
        assert!(comp.is_close(dzdx, 0.1));
        assert!(comp.is_close(dzdy, 0.1));
        assert!(comp.is_close(dzdx2, 0.0));
        assert!(comp.is_close(dzdy2, 0.0));
        assert!(comp.is_close(dzdxy, 0.0));

        let ze = spline2d
            .eval_extrap(4.0, 4.0, &mut xacc, &mut yacc)
            .unwrap();
        assert!(comp.is_close(ze, 1.8));
    }

    #[test]
    fn test_dyn_spline2d() {
        let xa = [0.0, 1.0, 2.0, 3.0];
        let ya = [0.0, 1.0, 2.0, 3.0];
        #[rustfmt::skip]
        let za = [
            1.0, 1.1, 1.2, 1.3,
            1.1, 1.2, 1.3, 1.4,
            1.2, 1.3, 1.4, 1.5,
            1.3, 1.4, 1.5, 1.6,
        ];

        Spline2d::new(Bicubic, &xa, &ya, &za).unwrap();
        Spline2d::new_dyn(Bicubic, &xa, &ya, &za).unwrap();
    }
}
