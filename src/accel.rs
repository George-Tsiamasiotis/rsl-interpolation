/// Index Look-up Acceration
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

    #[allow(dead_code)]
    #[doc(alias = "gsl_interp_bsearch")]
    /// This function returns the index i of the array `xarray` such that
    /// `xarray[i] <= x <= xarray[i+1]`. The index is searched for in the range [ilo, ihi].
    pub(crate) fn bsearch<T>(&self, xarray: &[T], x: T, ilo: usize, ihi: usize) -> usize
    where
        T: num::Float,
    {
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

    #[allow(dead_code)]
    #[doc(alias = "gsl_interp_accel_find")]
    /// Performs a lookup action on the data array. Returns an index i such that
    /// xarray[i] <= x < xarray[i+1].
    pub(crate) fn find<T>(&mut self, xarray: &[T], x: T) -> usize
    where
        T: num::Float,
    {
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
}

impl std::fmt::Debug for Accelerator {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Accelerator")
            .field("cache", &self.cache)
            .field("hits", &self.hits)
            .field("misses", &self.misses)
            .finish()
    }
}

impl Default for Accelerator {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_instantiation() {
        let _ = Accelerator::new();
    }

    #[test]
    fn test_reset() {
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
    fn test_debug_trait() {
        let acc = Accelerator::new();
        let _ = format!("{acc:?}");
        let _ = format!("{acc:#?}");
    }

    #[test]
    fn test_default_trait() {
        let _ = Accelerator::default();
    }
}
