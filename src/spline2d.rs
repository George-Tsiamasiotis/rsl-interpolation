//! Definition of the higher level `Spline2d` object.

use crate::{
    Accelerator2d, BuildInterpolator2d, Domain2dError, Interpolation2d, InterpolationError,
};

/// 2D Higher level interface.
///
/// A 2D Spline owns the data it is constructed with, and provides the same evaluation methods as the
/// lower-level Interpolator object, without needing to provide the data arrays in every call.
#[doc(alias = "gsl_spline2d")]
#[derive(Clone)]
pub struct Spline2d {
    /// The lower-level [`Interpolation2d`] object.
    interpolator: Box<dyn Interpolation2dClone>,
    /// The owned x data.
    xa: Box<[f64]>,
    /// The owned y data.
    ya: Box<[f64]>,
    /// The owned z data.
    za: Box<[f64]>,
}

impl Spline2d {
    /// Constructs a `Spline2d` of type `T` from the `xa`, `ya` and `za` arrays.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
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
    /// let interp = BilinearInterpolator::build(&xa, &ya, &za)?;
    /// let spline = Spline2d::from_interpolator(interp, &xa, &ya, &za);
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an [`InterpolationError`] if `T::build` fails. See implementor's documentation for
    /// possible errors.
    #[doc(alias = "gsl_spline2d_init")]
    pub fn build<T: BuildInterpolator2d + Clone>(
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
    ) -> Result<Self, InterpolationError> {
        let interp = T::build(xa, ya, za)?;
        Ok(Self::from_interpolator(interp, xa, ya, za))
    }

    /// Constructs a `Spline2d` from an Interpolator and its `xa`, `ya` and `za` arrays.
    ///
    /// # Example
    /// ```
    /// # use rsl_interpolation::*;
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
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
    /// let spline = Spline2d::build::<BicubicInterpolator>(&xa, &ya, &za)?;
    /// # Ok(())
    /// # }
    /// ```
    #[doc(alias = "gsl_spline_init")]
    #[must_use]
    pub fn from_interpolator(
        interpolator: impl Interpolation2d + Clone,
        xa: &[f64],
        ya: &[f64],
        za: &[f64],
    ) -> Self {
        Self {
            xa: xa.to_owned().into_boxed_slice(),
            ya: ya.to_owned().into_boxed_slice(),
            za: za.to_owned().into_boxed_slice(),
            interpolator: Box::new(interpolator),
        }
    }

    /// Returns a reference to the `Spline2d`'s `xa` data.
    #[must_use]
    pub fn xa(&self) -> &[f64] {
        &self.xa
    }

    /// Returns a reference to the `Spline2d`'s `ya` data.
    #[must_use]
    pub fn ya(&self) -> &[f64] {
        &self.ya
    }

    /// Returns a reference to the `Spline2d`'s `za` data.
    #[must_use]
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
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
    /// let spline = Spline2d::build::<BilinearInterpolator>(&xa, &ya, &za)?;
    ///
    /// let z = spline.eval(1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(z, 4.5);
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
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
    /// let spline = Spline2d::build::<BilinearInterpolator>(&xa, &ya, &za)?;
    ///
    /// let z = spline.eval_extrap(3.0, 6.0, &mut acc);
    /// assert_relative_eq!(z, 9.0);
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
    pub fn eval_extrap(&self, x: f64, y: f64, acc: &mut Accelerator2d) -> f64 {
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
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
    /// let spline = Spline2d::build::<BilinearInterpolator>(&xa, &ya, &za)?;
    ///
    /// let dzdx = spline.eval_deriv_x(1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(dzdx, 3.0);
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
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
    /// let spline = Spline2d::build::<BilinearInterpolator>(&xa, &ya, &za)?;
    ///
    /// let dzdx = spline.eval_deriv_y(1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(dzdx, 6.0);
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
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
    /// let spline = Spline2d::build::<BilinearInterpolator>(&xa, &ya, &za)?;
    ///
    /// let dzdx2 = spline.eval_deriv_xx(1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(dzdx2, 0.0); // Linear interpolation
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
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
    /// let spline = Spline2d::build::<BilinearInterpolator>(&xa, &ya, &za)?;
    ///
    /// let dzdx2 = spline.eval_deriv_yy(1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(dzdx2, 0.0); // Linear interpolation
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
    /// # use approx::assert_relative_eq;
    /// # fn main() -> Result<(), InterpolationError>{
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
    /// let spline = Spline2d::build::<BilinearInterpolator>(&xa, &ya, &za)?;
    ///
    /// let dzdxy = spline.eval_deriv_xy(1.5, 3.0, &mut acc)?;
    /// assert_relative_eq!(dzdxy, 0.0);
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

/// [`Box<dyn Interpolation2d>`] with Clone.
///
/// HACK: from
/// <https://stackoverflow.com/questions/30353462/how-to-clone-a-struct-storing-a-boxed-trait-object>.
trait Interpolation2dClone: Interpolation2d + Send + Sync {
    fn clone_box(&self) -> Box<dyn Interpolation2dClone>;
}

impl<T> Interpolation2dClone for T
where
    T: Clone + Interpolation2d,
{
    fn clone_box(&self) -> Box<dyn Interpolation2dClone> {
        Box::new(self.clone())
    }
}

impl Clone for Box<dyn Interpolation2dClone> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::*;

    #[test]
    fn build() {
        let xa = [0.0, 1.0, 2.0, 3.0];
        let ya = [0.0, 2.0, 4.0, 6.0];
        #[rustfmt::skip]
        let za = [
            0.0, 1.0, 2.0, 3.0,
            2.0, 3.0, 4.0, 5.0,
            4.0, 5.0, 6.0, 7.0,
            6.0, 7.0, 8.0, 9.0,
        ];

        let spline = Spline2d::build::<BicubicInterpolator>(&xa, &ya, &za).unwrap();
        let spline = spline.clone();

        assert_eq!(spline.xa(), xa);
        assert_eq!(spline.ya(), ya);
    }

    #[test]
    fn from_interpolator() {
        let xa = [0.0, 1.0, 2.0, 3.0];
        let ya = [0.0, 2.0, 4.0, 6.0];
        #[rustfmt::skip]
        let za = [
            0.0, 1.0, 2.0, 3.0,
            2.0, 3.0, 4.0, 5.0,
            4.0, 5.0, 6.0, 7.0,
            6.0, 7.0, 8.0, 9.0,
        ];

        let interp = BicubicInterpolator::build(&xa, &ya, &za).unwrap();
        let spline = Spline2d::from_interpolator(interp, &xa, &ya, &za);
        let spline = spline.clone();

        assert_eq!(spline.xa(), xa);
        assert_eq!(spline.ya(), ya);
    }
}
