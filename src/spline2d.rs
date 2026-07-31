//! Definition of the higher level Spline2d object.

use crate::Accelerator2d;
use crate::Domain2dError;
use crate::{DynInterpolator2d, Interpolation2d, Interpolation2dType};

/// 2D Higher level interface.
///
/// A 2D Spline owns the data it is constructed with, and provides the same evaluation methods as the
/// lower-level Interpolator object, without needing to provide the data arrays in every call.
#[doc(alias = "gsl_spline2d")]
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct Spline2d {
    /// The lower-level [`Interpolation2d`] object.
    interpolator: DynInterpolator2d,
    /// The owned x data.
    xa: Box<[f64]>,
    /// The owned y data.
    ya: Box<[f64]>,
    /// The owned z data.
    za: Box<[f64]>,
    /// The [`DynInterpolator2d`]'s [`Interpolation2dType`].
    typ: Interpolation2dType,
}

impl Spline2d {
    /// Constructs a 2D Spline from a [`DynInterpolator2d`] and its  `xa`, `ya` and `za` arrays.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let mut acc = Accelerator2d::new();
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
    /// let static_interp = BicubicInterpolator::build(&xa, &ya, &za)?;
    ///
    /// let dyn_interp = DynInterpolator2d::build(Interpolation2dType::Bicubic, &xa, &ya, &za)?;
    /// let spline = Spline2d::new(dyn_interp, &xa, &ya, &za);
    ///
    /// let x = 1.5;
    /// let y = 3.0;
    /// let y_interp = static_interp.eval(&xa, &ya, &za, x, y, &mut acc)?;
    /// let y_spline = spline.eval(x, y, &mut acc)?;
    ///
    /// assert_eq!(y_interp, y_spline);
    /// #
    /// # Ok(())
    /// # }
    /// ```
    #[doc(alias = "gsl_spline2d_init")]
    pub fn new(interpolator: DynInterpolator2d, xa: &[f64], ya: &[f64], za: &[f64]) -> Self {
        // in case the wrong arrays are passed
        if xa.len() * ya.len() != za.len() {
            panic!("Data array mismatch when constructing Spline")
        };
        Self {
            typ: interpolator.typ(),
            xa: xa.into(),
            ya: ya.into(),
            za: za.into(),
            interpolator,
        }
    }

    /// Returns the spline's [`Interpolation2dType`].
    pub fn typ(&self) -> Interpolation2dType {
        self.typ
    }

    /// Returns a reference to the spline's `xa` data.
    pub fn xa(&self) -> &[f64] {
        &self.xa
    }

    /// Returns a reference to the spline's `ya` data.
    pub fn ya(&self) -> &[f64] {
        &self.ya
    }

    /// Returns a reference to the spline's `za` data.
    pub fn za(&self) -> &[f64] {
        &self.za
    }

    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the
    /// [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let mut acc = Accelerator2d::new();
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
    /// let interp = DynInterpolator2d::build(Interpolation2dType::Bilinear, &xa, &ya, &za)?;
    /// let spline = Spline2d::new(interp, &xa, &ya, &za);
    /// #
    /// let z = spline.eval(1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(z, 4.5);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_spline2d_eval")]
    #[doc(alias = "gsl_spline2d_eval_e")]
    pub fn eval(&self, x: f64, y: f64, acc: &mut Accelerator2d) -> Result<f64, Domain2dError> {
        self.interpolator
            .eval(&self.xa, &self.ya, &self.za, x, y, acc)
    }

    /// Returns the interpolated value of `z` for a given point (`x`, `y`), using the
    /// [`Accelerator2d`] `acc`.
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
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let mut acc = Accelerator2d::new();
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
    /// let interp = DynInterpolator2d::build(Interpolation2dType::Bilinear, &xa, &ya, &za)?;
    /// let spline = Spline2d::new(interp, &xa, &ya, &za);
    /// #
    /// let z = spline.eval_extrap(3.0, 6.0, &mut acc)?;
    ///
    /// assert_eq!(z, 9.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_spline2d_eval_extrap")]
    #[doc(alias = "gsl_spline2d_eval_extrap_e")]
    pub fn eval_extrap(
        &self,
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator
            .eval_extrap(&self.xa, &self.ya, &self.za, x, y, acc)
    }

    /// Returns the interpolated value `d = ∂z/∂x` for a given point (`x`, `y`), using the
    /// [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let mut acc = Accelerator2d::new();
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
    /// let interp = DynInterpolator2d::build(Interpolation2dType::Bilinear, &xa, &ya, &za)?;
    /// let spline = Spline2d::new(interp, &xa, &ya, &za);
    ///
    /// let dzdx = spline.eval_deriv_x(1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(dzdx, 3.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_spline2d_eval_deriv_x")]
    #[doc(alias = "gsl_spline2d_eval_deriv_x_e")]
    pub fn eval_deriv_x(
        &self,
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator
            .eval_deriv_x(&self.xa, &self.ya, &self.za, x, y, acc)
    }

    /// Returns the interpolated value `d = ∂z/∂y` for a given point (`x`, `y`), using the
    /// [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let mut acc = Accelerator2d::new();
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
    /// let interp = DynInterpolator2d::build(Interpolation2dType::Bilinear, &xa, &ya, &za)?;
    /// let spline = Spline2d::new(interp, &xa, &ya, &za);
    ///
    /// let dzdx = spline.eval_deriv_y(1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(dzdx, 6.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_spline2d_eval_deriv_y")]
    #[doc(alias = "gsl_spline2d_eval_deriv_y_e")]
    pub fn eval_deriv_y(
        &self,
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator
            .eval_deriv_y(&self.xa, &self.ya, &self.za, x, y, acc)
    }

    /// Returns the interpolated value `d = 𝜕²z/𝜕x²` for a given point (`x`, `y`), using the
    /// [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let mut acc = Accelerator2d::new();
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
    /// let interp = DynInterpolator2d::build(Interpolation2dType::Bilinear, &xa, &ya, &za)?;
    /// let spline = Spline2d::new(interp, &xa, &ya, &za);
    ///
    /// let dzdx2 = spline.eval_deriv_xx(1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(dzdx2, 0.0); // Linear interpolation
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_interp2d_eval_deriv_xx")]
    #[doc(alias = "gsl_interp2d_eval_deriv_xx_e")]
    pub fn eval_deriv_xx(
        &self,
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator
            .eval_deriv_xx(&self.xa, &self.ya, &self.za, x, y, acc)
    }

    /// Returns the interpolated value `d = 𝜕²z/𝜕y²` for a given point (`x`, `y`), using the
    /// [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let mut acc = Accelerator2d::new();
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
    ///
    /// let interp = DynInterpolator2d::build(Interpolation2dType::Bilinear, &xa, &ya, &za)?;
    /// let spline = Spline2d::new(interp, &xa, &ya, &za);
    ///
    /// let dzdx2 = spline.eval_deriv_yy(1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(dzdx2, 0.0); // Linear interpolation
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_interp2d_eval_deriv_yy")]
    #[doc(alias = "gsl_interp2d_eval_deriv_yy_e")]
    pub fn eval_deriv_yy(
        &self,
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator
            .eval_deriv_yy(&self.xa, &self.ya, &self.za, x, y, acc)
    }

    /// Returns the interpolated value `d = 𝜕²z/𝜕x𝜕y` for a given point (`x`, `y`), using the
    /// [`Accelerator2d`] `acc`.
    ///
    /// # Example
    ///
    /// ```
    /// # use rsl_interpolation::*;
    /// #
    /// # fn main() -> Result<(), InterpolatorError>{
    /// let mut acc = Accelerator2d::new();
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
    /// let interp = DynInterpolator2d::build(Interpolation2dType::Bilinear, &xa, &ya, &za)?;
    /// let spline = Spline2d::new(interp, &xa, &ya, &za);
    ///
    /// let dzdxy = spline.eval_deriv_xy(1.5, 3.0, &mut acc)?;
    ///
    /// assert_eq!(dzdxy, 0.0);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns a [`Domain2dError`] if `x` is outside the range of `xa` or `y` is outside the range
    /// of `ya`.
    #[doc(alias = "gsl_interp2d_eval_deriv_xy")]
    #[doc(alias = "gsl_interp2d_eval_deriv_xy_e")]
    pub fn eval_deriv_xy(
        &self,
        x: f64,
        y: f64,
        acc: &mut Accelerator2d,
    ) -> Result<f64, Domain2dError> {
        self.interpolator
            .eval_deriv_xy(&self.xa, &self.ya, &self.za, x, y, acc)
    }
}
