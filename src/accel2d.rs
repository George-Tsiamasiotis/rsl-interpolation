//! 2D Index Look-up Acceleration object.

use crate::Accelerator;
use crate::z_idx;

/// 2D Index Look-up Acceleration.
///
/// This object caches the values extracted from the `za` array, which are more likely to be used
/// again, either by evaluating the spline to a nearby point, or/and by evaluating it's derivates
/// as well.
///
/// Only the **grid** points are cached, **not any interpolated values**.
///
/// This optimization seems to be quite effective (up to +50% from my testing), and could be
/// important in cases where `za` arrays are big and don't fit in the CPU's cache. It is of course,
/// very situational, since it also depends on the way one evaluates his splines.
///
/// The overhead of cache misses should be truly negligible, since the process just falls back to
/// calculating the values in the usual manner.
#[derive(Debug, Clone)]
pub struct Accelerator2d {
    xacc: Accelerator,
    yacc: Accelerator,
    xi: usize,
    yi: usize,
    xgrid_values: (f64, f64),
    ygrid_values: (f64, f64),
    zgrid_values: (f64, f64, f64, f64),
    partials: (f64, f64),
}

impl Accelerator2d {
    /// Creates a new empty `Accelerator2d`.
    #[must_use]
    pub fn new() -> Self {
        Self {
            xacc: Accelerator::new(),
            yacc: Accelerator::new(),
            xi: usize::MAX,
            yi: usize::MAX,
            xgrid_values: (f64::NAN, f64::NAN),
            ygrid_values: (f64::NAN, f64::NAN),
            zgrid_values: (f64::NAN, f64::NAN, f64::NAN, f64::NAN),
            partials: (f64::NAN, f64::NAN),
        }
    }

    /// Resets the `Accelerator2d`.
    pub fn reset(&mut self) {
        *self = Self::new()
    }

    /// Returns a mutable reference to the x-data 1D [`Accelerator`].
    pub fn xacc(&mut self) -> &mut Accelerator {
        &mut self.xacc
    }

    /// Returns a mutable reference to the y-data 1D [`Accelerator`].
    pub fn yacc(&mut self) -> &mut Accelerator {
        &mut self.yacc
    }
}

/// Update the cache's fields. Should be called in a specific order.
impl Accelerator2d {
    pub(crate) fn update(&mut self, xa: &[f64], ya: &[f64], za: &[f64], x: f64, y: f64) {
        let old_xi = self.xi;
        let old_yi = self.yi;
        self.xi = self.xacc.find(xa, x);
        self.yi = self.yacc.find(ya, y);

        // Update grid values only if acc indices changed
        if (self.xi == old_xi) && (self.yi == old_yi) {
            // Make sure that grid values are initialized. Uninitialized fields may appear when 1D
            // interpolators use the Accelerator before a 2D interpolator does.
            #[expect(clippy::redundant_else, reason = "clearer")]
            if self.zgrid_values.0.is_finite() {
                return;
            } else {
                std::hint::cold_path();
            }
        }
        let xi = self.xi;
        let yi = self.yi;
        let xlen = xa.len();
        let ylen = ya.len();
        self.xgrid_values = (xa[xi], xa[xi + 1]);
        self.ygrid_values = (ya[yi], ya[yi + 1]);
        self.zgrid_values.0 = za[z_idx(xi, yi, xlen, ylen)];
        self.zgrid_values.1 = za[z_idx(xi, yi + 1, xlen, ylen)];
        self.zgrid_values.2 = za[z_idx(xi + 1, yi, xlen, ylen)];
        self.zgrid_values.3 = za[z_idx(xi + 1, yi + 1, xlen, ylen)];
        self.partials.0 = self.xgrid_values.1 - self.xgrid_values.0;
        self.partials.1 = self.ygrid_values.1 - self.ygrid_values.0;
    }
}

/// Getter methods for grid point quantities.
impl Accelerator2d {
    pub(crate) fn get_xy_grid_values(&self) -> (f64, f64, f64, f64) {
        (
            self.xgrid_values.0,
            self.xgrid_values.1,
            self.ygrid_values.0,
            self.ygrid_values.1,
        )
    }

    pub(crate) fn get_z_grid_values(&self) -> (f64, f64, f64, f64) {
        self.zgrid_values
    }

    pub(crate) fn get_partials(&self) -> (f64, f64) {
        self.partials
    }
}

impl Default for Accelerator2d {
    fn default() -> Self {
        Self::new()
    }
}
