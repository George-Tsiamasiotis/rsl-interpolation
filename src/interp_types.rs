/// Interpolation types supported by GSL.
///
/// ## Example
/// ```
/// # use rsl_interpolation::InterpolationType;
/// #
/// # fn main() {
/// let interp_type = InterpolationType::Akima;
/// # }
/// ```
/// Descriptions are taken directly from GSL's
/// [interpolation](https://www.gnu.org/software/gsl/doc/html/interp.html#c.gsl_interp_type) page
#[derive(Debug)]
pub enum InterpolationType {
    /// Linear interpolation. This interpolation method does not require any additional memory.
    Linear,

    /// Polynomial interpolation. This method should only be used for interpolating small number
    /// of points because polynomial interpolation introduces large oscillations, even for
    /// well-behaved datasets. The number of terms in the interpolating polynomial is equal to
    /// the number of points.
    Polynomial,

    /// Cubic spline with natural boundary conditions. The resulting curve is piecewise cubic on
    /// each interval, with matching first and second derivatives at the supplied data-points.
    /// The second derivative is chosen to be zero at the first point and last point.
    Cubic,

    /// Cubic spline with periodic boundary conditions. The resulting curve is piecewise cubic
    /// on each interval, with matching first and second derivatives at the supplied data-points.
    /// The derivatives at the first and last points are also matched. Note that the last point
    /// in the data must have the same y-value as the first point, otherwise the resulting
    /// periodic interpolation will have a discontinuity at the boundary.
    CubicPeriodic,

    /// Non-rounded Akima spline with natural boundary conditions. This method uses the
    /// non-rounded corner algorithm of Wodicka.
    Akima,

    /// Non-rounded Akima spline with periodic boundary conditions. This method uses the
    /// non-rounded corner algorithm of Wodicka.
    AkimaPeriodic,
    ///Steffen’s method guarantees the monotonicity of the interpolating function between the
    ///given data points. Therefore, minima and maxima can only occur exactly at the data points,
    ///and there can never be spurious oscillations between data points. The interpolation
    /// function is piecewise cubic in each interval. The resulting curve and its first derivative
    /// are guaranteed to be continuous, but the second derivative may be discontinuous.
    Steffen,
}
