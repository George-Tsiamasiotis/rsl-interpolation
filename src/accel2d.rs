//! 2D Index Look-up Acceleration object.

use crate::Accelerator;
use crate::z_idx;

#[derive(Debug, Clone, Copy)]
/// 2D Index Look-up Acceleration
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
pub struct Accelerator2d {
    xacc: Accelerator,
    yacc: Accelerator,
    xgrid_values: (f64, f64),
    ygrid_values: (f64, f64),
    zgrid_values: (f64, f64, f64, f64),
    zx_values: (f64, f64, f64, f64),
    zy_values: (f64, f64, f64, f64),
    zxy_values: (f64, f64, f64, f64),
    partials: (f64, f64),
    uninit: bool,
}

impl Accelerator2d {
    /// Creates a new empty [`Accelerator2d`]
    pub fn new() -> Self {
        Self {
            xacc: Accelerator::new(),
            yacc: Accelerator::new(),
            xgrid_values: (f64::NAN, f64::NAN),
            ygrid_values: (f64::NAN, f64::NAN),
            zgrid_values: (f64::NAN, f64::NAN, f64::NAN, f64::NAN),
            zx_values: (f64::NAN, f64::NAN, f64::NAN, f64::NAN),
            zy_values: (f64::NAN, f64::NAN, f64::NAN, f64::NAN),
            zxy_values: (f64::NAN, f64::NAN, f64::NAN, f64::NAN),
            partials: (f64::NAN, f64::NAN),
            uninit: true,
        }
    }

    /// Resets the Accelerator.
    pub fn reset(&mut self) {
        *self = Self::new()
    }

    /// Resets the indices. Useful for benchmarking, to avoid the overhead of resetting all the
    /// fields at each iteration.
    pub fn soft_reset(&mut self) {
        self.xacc.cache = 0;
        self.yacc.cache = 0;
    }

    /// Returns a mutable reference to the x-data 1D [`Accelerator`].
    pub fn xacc(&mut self) -> &mut Accelerator {
        &mut self.xacc
    }

    /// Returns a mutable reference to the y-data 1D [`Accelerator`].
    pub fn yacc(&mut self) -> &mut Accelerator {
        &mut self.yacc
    }

    /// Checks if the Accelerator's cached values are up to date.
    pub(crate) fn is_uptodate(&mut self, xa: &[f64], ya: &[f64], x: f64, y: f64) -> bool {
        // The first time that the Cache is being called, the values are uninitialized, but the
        // interpolator does not know that. This forces the Cache to be updated the first time it
        // is called after initialization.
        //
        // It is also necessary to setup the indices and the Accelerators for the next evaluation.
        //
        // Every evaluation after that is not affected.

        let old_xi = self.xacc.cache;
        let old_yi = self.yacc.cache;
        self.xacc.find(xa, x);
        self.yacc.find(ya, y);

        if self.uninit {
            self.uninit = false;
            return false;
        }

        // If xa landed on the same interval (cache hit), `find()` simply returns the the cached
        // index, otherwise it recalculates it.

        (old_xi == self.xacc.cache) && (old_yi == self.yacc.cache)
    }

    pub(crate) fn update_step1(&mut self, xa: &[f64], ya: &[f64], za: &[f64]) {
        self.update_xy_grid_values(xa, ya);
        self.update_z_grid_values(za, xa.len(), ya.len());
        self.update_partials();
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn update_step2(
        &mut self,
        xa: &[f64],
        ya: &[f64],
        zx: &[f64],
        zy: &[f64],
        zxy: &[f64],
        dt: f64,
        du: f64,
    ) {
        self.update_zxminmaxing(zx, xa.len(), ya.len(), dt);
        self.update_zyminmaxing(zy, xa.len(), ya.len(), du);
        self.update_zxyminmaxing(zxy, xa.len(), ya.len(), dt * du);
    }
}

/// Methods for updating the cache's fields. Should be called in a specific order.
impl Accelerator2d {
    fn update_xy_grid_values(&mut self, xa: &[f64], ya: &[f64]) {
        let xi = self.xacc.cache;
        let yi = self.yacc.cache;
        self.xgrid_values = (xa[xi], xa[xi + 1]);
        self.ygrid_values = (ya[yi], ya[yi + 1]);
    }

    fn update_z_grid_values(&mut self, za: &[f64], xlen: usize, ylen: usize) {
        let xi = self.xacc.cache;
        let yi = self.yacc.cache;
        self.zgrid_values.0 = za[z_idx(xi, yi, xlen, ylen)];
        self.zgrid_values.1 = za[z_idx(xi, yi + 1, xlen, ylen)];
        self.zgrid_values.2 = za[z_idx(xi + 1, yi, xlen, ylen)];
        self.zgrid_values.3 = za[z_idx(xi + 1, yi + 1, xlen, ylen)];
    }

    fn update_partials(&mut self) {
        self.partials.0 = self.xgrid_values.1 - self.xgrid_values.0;
        self.partials.1 = self.ygrid_values.1 - self.ygrid_values.0;
    }

    fn update_zxminmaxing(&mut self, zx: &[f64], xsize: usize, ysize: usize, d: f64) {
        let xi = self.xacc.cache;
        let yi = self.yacc.cache;
        self.zx_values.0 = zx[z_idx(xi, yi, xsize, ysize)] / d;
        self.zx_values.1 = zx[z_idx(xi, yi + 1, xsize, ysize)] / d;
        self.zx_values.2 = zx[z_idx(xi + 1, yi, xsize, ysize)] / d;
        self.zx_values.3 = zx[z_idx(xi + 1, yi + 1, xsize, ysize)] / d;
    }

    fn update_zyminmaxing(&mut self, zy: &[f64], xsize: usize, ysize: usize, d: f64) {
        let xi = self.xacc.cache;
        let yi = self.yacc.cache;
        self.zy_values.0 = zy[z_idx(xi, yi, xsize, ysize)] / d;
        self.zy_values.1 = zy[z_idx(xi, yi + 1, xsize, ysize)] / d;
        self.zy_values.2 = zy[z_idx(xi + 1, yi, xsize, ysize)] / d;
        self.zy_values.3 = zy[z_idx(xi + 1, yi + 1, xsize, ysize)] / d;
    }

    fn update_zxyminmaxing(&mut self, zxy: &[f64], xsize: usize, ysize: usize, d: f64) {
        let xi = self.xacc.cache;
        let yi = self.yacc.cache;
        self.zxy_values.0 = zxy[z_idx(xi, yi, xsize, ysize)] / d;
        self.zxy_values.1 = zxy[z_idx(xi, yi + 1, xsize, ysize)] / d;
        self.zxy_values.2 = zxy[z_idx(xi + 1, yi, xsize, ysize)] / d;
        self.zxy_values.3 = zxy[z_idx(xi + 1, yi + 1, xsize, ysize)] / d;
    }
}

/// Getter methods for grid point quantities
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

    pub(crate) fn get_zxminmaxxing(&self) -> (f64, f64, f64, f64) {
        self.zx_values
    }

    pub(crate) fn get_zyminmaxxing(&self) -> (f64, f64, f64, f64) {
        self.zy_values
    }

    pub(crate) fn get_zxyminmaxxing(&self) -> (f64, f64, f64, f64) {
        self.zxy_values
    }
}

impl Default for Accelerator2d {
    fn default() -> Self {
        Self::new()
    }
}
