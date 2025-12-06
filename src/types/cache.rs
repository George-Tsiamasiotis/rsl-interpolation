use crate::Accelerator;
use crate::DomainError;
use crate::z_idx;

#[derive(Debug, Clone, Copy)]
/// Cache holding values retrieved from the `za` array.
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
/// The overhead of cache misses should be trully negligible, since the process just falls back to
/// calculating the values in the usual manner.
pub struct Cache<T> {
    acc_indices: (usize, usize),
    xgrid_values: (T, T),
    ygrid_values: (T, T),
    zgrid_values: (T, T, T, T),
    zx_values: (T, T, T, T),
    zy_values: (T, T, T, T),
    zxy_values: (T, T, T, T),
    partials: (T, T),
    uninit: bool,
}

impl<T> Cache<T>
where
    T: crate::Num,
{
    /// Creates a new empty [`Cache2d`]
    pub fn new() -> Self {
        let def = T::nan();
        Self {
            acc_indices: (0, 0),
            xgrid_values: (def, def),
            ygrid_values: (def, def),
            zgrid_values: (def, def, def, def),
            zx_values: (def, def, def, def),
            zy_values: (def, def, def, def),
            zxy_values: (def, def, def, def),
            partials: (def, def),
            uninit: true,
        }
    }

    pub fn reset(&mut self) {
        *self = Self::new()
    }

    /// Resets the indeces. Useful for benchmarking, to avoid the overhead of resetting all the
    /// fields at each iteration.
    pub fn soft_reset(&mut self) {
        self.acc_indices = (0, 0)
    }

    pub(crate) fn is_uptodate(&mut self, xa: &[T], ya: &[T], x: T, y: T) -> bool {
        // The first time that the Cache is being called, the values are uninitialized, but the
        // interpolator does not know that. This forces the Cache to be updated the first time it
        // is called after initialization.
        //
        // Every evaluation after that is not affected.
        if self.uninit {
            self.uninit = false;
            return false;
        }

        let xi = self.acc_indices.0;
        let yi = self.acc_indices.1;
        let x_inbounds: bool = (x > xa[xi]) && (x < xa[xi + 1]);
        let y_inbounds: bool = (y > ya[yi]) && (y < ya[yi + 1]);
        x_inbounds && y_inbounds
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn update_step1(
        &mut self,
        xa: &[T],
        ya: &[T],
        za: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) -> Result<(), DomainError> {
        self.update_acc_indeces(xa, ya, x, y, xacc, yacc);
        self.update_xy_grid_values(xa, ya);
        self.update_z_grid_values(za, xa.len(), ya.len())?;
        self.update_partials();
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn update_step2(
        &mut self,
        xa: &[T],
        ya: &[T],
        zx: &[T],
        zy: &[T],
        zxy: &[T],
        dt: T,
        du: T,
    ) -> Result<(), DomainError> {
        self.update_zxminmaxing(zx, xa.len(), ya.len(), dt)?;
        self.update_zyminmaxing(zy, xa.len(), ya.len(), du)?;
        self.update_zxyminmaxing(zxy, xa.len(), ya.len(), dt * du)?;
        Ok(())
    }
}

/// Methods for updating the cache's fields. Should be called in a specific order.
impl<T> Cache<T>
where
    T: crate::Num,
{
    fn update_acc_indeces(
        &mut self,
        xa: &[T],
        ya: &[T],
        x: T,
        y: T,
        xacc: &mut Accelerator,
        yacc: &mut Accelerator,
    ) {
        self.acc_indices.0 = xacc.find(xa, x);
        self.acc_indices.1 = yacc.find(ya, y);
    }

    fn update_xy_grid_values(&mut self, xa: &[T], ya: &[T]) {
        let xi = self.acc_indices.0;
        let yi = self.acc_indices.1;
        self.xgrid_values = (xa[xi], xa[xi + 1]);
        self.ygrid_values = (ya[yi], ya[yi + 1]);
    }

    fn update_z_grid_values(
        &mut self,
        za: &[T],
        xlen: usize,
        ylen: usize,
    ) -> Result<(), DomainError> {
        let xi = self.acc_indices.0;
        let yi = self.acc_indices.1;
        self.zgrid_values.0 = za[z_idx(xi, yi, xlen, ylen)?];
        self.zgrid_values.1 = za[z_idx(xi, yi + 1, xlen, ylen)?];
        self.zgrid_values.2 = za[z_idx(xi + 1, yi, xlen, ylen)?];
        self.zgrid_values.3 = za[z_idx(xi + 1, yi + 1, xlen, ylen)?];
        Ok(())
    }

    fn update_partials(&mut self) {
        self.partials.0 = self.xgrid_values.1 - self.xgrid_values.0;
        self.partials.1 = self.ygrid_values.1 - self.ygrid_values.0;
    }

    fn update_zxminmaxing(
        &mut self,
        zx: &[T],
        xsize: usize,
        ysize: usize,
        d: T,
    ) -> Result<(), DomainError> {
        let xi = self.acc_indices.0;
        let yi = self.acc_indices.1;
        self.zx_values.0 = zx[z_idx(xi, yi, xsize, ysize)?] / d;
        self.zx_values.1 = zx[z_idx(xi, yi + 1, xsize, ysize)?] / d;
        self.zx_values.2 = zx[z_idx(xi + 1, yi, xsize, ysize)?] / d;
        self.zx_values.3 = zx[z_idx(xi + 1, yi + 1, xsize, ysize)?] / d;
        Ok(())
    }

    fn update_zyminmaxing(
        &mut self,
        zy: &[T],
        xsize: usize,
        ysize: usize,
        d: T,
    ) -> Result<(), DomainError> {
        let xi = self.acc_indices.0;
        let yi = self.acc_indices.1;
        self.zy_values.0 = zy[z_idx(xi, yi, xsize, ysize)?] / d;
        self.zy_values.1 = zy[z_idx(xi, yi + 1, xsize, ysize)?] / d;
        self.zy_values.2 = zy[z_idx(xi + 1, yi, xsize, ysize)?] / d;
        self.zy_values.3 = zy[z_idx(xi + 1, yi + 1, xsize, ysize)?] / d;
        Ok(())
    }

    fn update_zxyminmaxing(
        &mut self,
        zxy: &[T],
        xsize: usize,
        ysize: usize,
        d: T,
    ) -> Result<(), DomainError> {
        let xi = self.acc_indices.0;
        let yi = self.acc_indices.1;
        self.zxy_values.0 = zxy[z_idx(xi, yi, xsize, ysize)?] / d;
        self.zxy_values.1 = zxy[z_idx(xi, yi + 1, xsize, ysize)?] / d;
        self.zxy_values.2 = zxy[z_idx(xi + 1, yi, xsize, ysize)?] / d;
        self.zxy_values.3 = zxy[z_idx(xi + 1, yi + 1, xsize, ysize)?] / d;
        Ok(())
    }
}

/// Getter methods for grid point quantities
impl<T> Cache<T>
where
    T: crate::Num,
{
    pub(crate) fn get_xy_indeces(&self) -> (usize, usize) {
        (self.acc_indices.0, self.acc_indices.1)
    }

    pub(crate) fn get_xy_grid_values(&self) -> (T, T, T, T) {
        (
            self.xgrid_values.0,
            self.xgrid_values.1,
            self.ygrid_values.0,
            self.ygrid_values.1,
        )
    }

    pub(crate) fn get_z_grid_values(&self) -> (T, T, T, T) {
        self.zgrid_values
    }

    pub(crate) fn get_partials(&self) -> (T, T) {
        self.partials
    }

    pub(crate) fn get_zxminmaxxing(&self) -> (T, T, T, T) {
        self.zx_values
    }

    pub(crate) fn get_zyminmaxxing(&self) -> (T, T, T, T) {
        self.zy_values
    }

    pub(crate) fn get_zxyminmaxxing(&self) -> (T, T, T, T) {
        self.zxy_values
    }
}

impl<T> Default for Cache<T>
where
    T: crate::Num,
{
    fn default() -> Self {
        Self::new()
    }
}
