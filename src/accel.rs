/// Index Look-up Acceration
///
/// This object caches the previous value of an index lookup. When the subsequent interpolation
/// point falls in the same interval, its index value can be returned immediately.
///
/// The performance boost can be significant when continuously evalulating splines around the same
/// area as the previous point. Moreover, the same Accelerator can be shared across multiple
/// Splines, if they are defined over the same x points. This is especially useful ODE solvers,
/// as the solver's step size is usually much smaller that the xarray spacing.
///
/// See [`GSL's Acceleration section`].
///
/// ## Example
///
/// ```
/// # use rsl_interpolation::Accelerator;
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
    /// This function returns the index i of the array `xarray` such that
    /// `xarray[i] <= x <= xarray[i+1]`. The index is searched for in the range [ilo, ihi].
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

    /// Resets the Accelerator's stats to 0.
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
        let _ = format!("{:?}", acc);
        let _ = format!("{:#?}", acc);
    }
}
