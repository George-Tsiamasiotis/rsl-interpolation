//! 1D Index Look-up Acceleration object.

/// 1D Index Look-up Acceleration.
///
/// This object caches the previous value of an index lookup. When the subsequent interpolation
/// point falls in the same interval, its index value can be returned immediately.
///
/// The performance boost can be significant when continuously evaluating splines around the same
/// area as the previous point. Moreover, the same Accelerator can be shared across multiple
/// Splines, if they are defined over the same x points. This is especially useful in ODE solvers,
/// as the solver's step size is usually much smaller that the xarray spacing.
///
/// See [`GSL's Acceleration section`].
///
/// ## Example
///
/// ```
/// # use rsl_interpolation::*;
/// #
/// # fn main() {
/// let mut acc = Accelerator::new();
/// # }
/// ```
/// [`GSL's Acceleration section`]: https://www.gnu.org/software/gsl/doc/html/interp.html#d-index-look-up-and-acceleration
#[doc(alias = "gsl_interp_accel")]
#[derive(Default, Debug, Clone, Copy)]
pub struct Accelerator {
    /// The current cached index.
    pub(crate) cache: usize,
    /// The total cache hits.
    pub(crate) hits: usize,
    /// The total cache misses.
    pub(crate) misses: usize,
}

impl Accelerator {
    /// Creates a new Accelerator.
    pub fn new() -> Self {
        Accelerator {
            cache: 0,
            hits: 0,
            misses: 0,
        }
    }

    /// This function returns the index i of the array `xarray` such that
    /// `xarray[i] <= x <= xarray[i+1]`. The index is searched for in the range [ilo, ihi].
    #[doc(alias = "gsl_interp_bsearch")]
    pub(crate) fn bsearch(&self, xarray: &[f64], x: f64, ilo: usize, ihi: usize) -> usize {
        let mut ilo = ilo;
        let mut ihi = ihi;
        while ihi > ilo + 1 {
            let i = (ihi + ilo) / 2;
            if xarray[i] > x {
                ihi = i;
            } else {
                ilo = i;
            }
        }
        ilo
    }

    /// Performs a lookup action on the data array. Returns an index i such that
    /// xarray[i] <= x < xarray[i+1].
    #[doc(alias = "gsl_interp_accel_find")]
    pub(crate) fn find(&mut self, xarray: &[f64], x: f64) -> usize {
        if x < xarray[self.cache] {
            self.misses += 1;
            self.cache = self.bsearch(xarray, x, 0, self.cache);
        } else if x >= xarray[self.cache + 1] {
            self.misses += 1;
            self.cache = self.bsearch(xarray, x, self.cache, xarray.len() - 1);
        } else {
            self.hits += 1;
        }
        self.cache
    }

    /// Resets the Accelerator's stats to 0.
    #[doc(alias = "gsl_interp_accel_reset")]
    pub fn reset(&mut self) {
        self.cache = 0;
        self.hits = 0;
        self.misses = 0;
    }

    /// Returns the total cache hits of the Accelerator.
    pub fn hits(&self) -> usize {
        self.hits
    }

    /// Returns the total cache misses of the Accelerator.
    pub fn misses(&self) -> usize {
        self.misses
    }

    /// Returns the cached index of the Accelerator.
    pub fn cache(&self) -> usize {
        self.cache
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn setup_acc() -> Accelerator {
        Accelerator::new()
    }

    fn setup_xarray() -> [f64; 5] {
        [0.0, 1.0, 2.0, 3.0, 4.0]
    }

    #[test]
    fn instantiation() {
        let acc = Accelerator::default();
        let _ = format!("{:?}", acc);
        assert_eq!(acc.hits(), 0);
        assert_eq!(acc.misses(), 0);
        assert_eq!(acc.cache(), 0);
    }

    #[test]
    fn reset() {
        let mut acc = Accelerator::new();
        acc.cache = 1;
        acc.hits = 1;
        acc.misses = 1;
        acc.reset();

        assert_eq!(acc.cache, 0);
        assert_eq!(acc.hits, 0);
        assert_eq!(acc.misses, 0);
    }

    #[test]
    fn bsearch_interior_point() {
        let (acc, xarray) = (setup_acc(), setup_xarray());

        let res = acc.bsearch(&xarray, 1.5, 0, 4);
        assert_eq!(res, 1);
    }

    #[test]
    fn search_last_value() {
        let (acc, xarray) = (setup_acc(), setup_xarray());

        let res = acc.bsearch(&xarray, 4.0, 0, 4);
        assert_eq!(res, 3);
    }

    #[test]
    fn search_first_value() {
        let (acc, xarray) = (setup_acc(), setup_xarray());

        let res = acc.bsearch(&xarray, 0.0, 0, 4);
        assert_eq!(res, 0);
    }

    #[test]
    fn bsearch_boundary() {
        let (acc, xarray) = (setup_acc(), setup_xarray());

        let res = acc.bsearch(&xarray, 2.0, 0, 4);
        assert_eq!(res, 2);
    }

    #[test]
    fn bsearch_above_bounds() {
        let (acc, xarray) = (setup_acc(), setup_xarray());

        let res = acc.bsearch(&xarray, 10.0, 0, 4);
        assert_eq!(res, 3);
    }

    #[test]
    fn bsearch_below_bounds() {
        let (acc, xarray) = (setup_acc(), setup_xarray());

        let res = acc.bsearch(&xarray, -10.0, 0, 4);
        assert_eq!(res, 0);
    }

    #[test]
    fn accelerator() {
        let xarray = setup_xarray();
        let mut acc = setup_acc();
        let mut k1 = 0;
        let mut k2 = 0;
        let mut t = false;
        let r = [
            -0.2, 0.0, 0.1, 0.7, 1.0, 1.3, 1.9, 2.0, 2.2, 2.7, 3.0, 3.1, 3.6, 4.0, 4.1, 4.9,
        ];

        // run through all the pairs of points.
        while (k1 < 16) & (k2 < 16) {
            let x = if t { r[k1] } else { r[k2] };
            t = !t;

            if !t {
                k1 = (k1 + 1) % 16;
                if k1 == 0 {
                    k2 += 1;
                }
            }

            let i = acc.find(&xarray, x);
            let j = acc.bsearch(&xarray, x, 0, 4);
            assert_eq!(i, j);
        }
    }
}
