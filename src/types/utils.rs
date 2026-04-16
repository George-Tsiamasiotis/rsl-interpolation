//! Common utility methods for the Interpolators.

use crate::{Domain1dError, Domain2dError, InterpolationError};

/// Checks that supplied datasets are valid.
pub(crate) fn check1d_data(
    xa: &[f64],
    ya: &[f64],
    min_size: usize,
) -> Result<(), InterpolationError> {
    if !xa.iter().is_sorted() {
        return Err(InterpolationError::UnsortedDataset);
    }
    if xa.len() != ya.len() {
        return Err(InterpolationError::DatasetMismatch);
    }
    if xa.len() < min_size {
        return Err(InterpolationError::NotEnoughPoints);
    }
    Ok(())
}

/// Checks that the passed xa, ya and za datasets are valid.
pub(crate) fn check2d_data(
    xa: &[f64],
    ya: &[f64],
    za: &[f64],
    min_size: usize,
) -> Result<(), InterpolationError> {
    if xa.len() * ya.len() != za.len() {
        return Err(InterpolationError::ZGridMismatch);
    }
    if (xa.len() < min_size) | (ya.len() < min_size) {
        return Err(InterpolationError::NotEnoughPoints);
    }
    if (!xa.iter().is_sorted()) | (!ya.iter().is_sorted()) {
        Err(InterpolationError::UnsortedDataset)
    } else {
        Ok(())
    }
}

pub(crate) fn check_if_inbounds(xa: &[f64], x: f64) -> Result<(), Domain1dError> {
    match (xa.first(), xa.last()) {
        (Some(first), Some(last)) if x >= *first && x <= *last => Ok(()),
        _ => Err(Domain1dError { x }),
    }
}

pub(crate) fn check_if_inbounds2d(
    xa: &[f64],
    ya: &[f64],
    x: f64,
    y: f64,
) -> Result<(), Domain2dError> {
    match (check_if_inbounds(xa, x), check_if_inbounds(ya, y)) {
        (Ok(_), Ok(_)) => Ok(()),
        _ => Err(Domain2dError { x, y }),
    }
}

/// Calculates the n-th discrete difference: out[i] = s[i+1]-s[i].
pub(crate) fn diff(s: &[f64]) -> Vec<f64> {
    s.windows(2)
        .map(|xy| {
            let [x, y] = xy else { unreachable!() };
            *y - *x
        })
        .collect()
}

/// Function for doing the spline integral evaluation, which is common to both the cspline and
/// akima methods.
pub(crate) fn integ_eval(ai: f64, bi: f64, ci: f64, di: f64, xi: f64, a: f64, b: f64) -> f64 {
    let r1 = a - xi;
    let r2 = b - xi;
    let r12 = r1.powi(2);
    let r22 = r2.powi(2);
    let rsum = r1 + r2;
    let bterm = 0.5 * bi * rsum;
    let cterm = (1.0 / 3.0) * ci * (r12 + r22 + r1 * r2);
    let dterm = 0.25 * di * rsum * (r12 + r22);

    (b - a) * (ai + bterm + cterm + dterm)
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn check_data() {
        let ya = [0.0, 1.0, 2.0];

        let xa = [0.0, 1.0, 2.0];
        assert!(check1d_data(&xa, &ya, 2).is_ok());
        assert!(matches!(
            check1d_data(&xa, &ya, 4).unwrap_err(),
            InterpolationError::NotEnoughPoints
        ));

        let xa = [2.0, 1.0, 2.0];
        assert!(matches!(
            check1d_data(&xa, &ya, 2).unwrap_err(),
            InterpolationError::UnsortedDataset
        ));

        let xa = [0.0, 1.0, 2.0, 3.0];
        assert!(matches!(
            check1d_data(&xa, &ya, 2).unwrap_err(),
            InterpolationError::DatasetMismatch
        ));
    }
    #[test]
    fn inbounds() {
        let xa = [0.0, 1.0, 2.0];

        assert!(check_if_inbounds(&xa, 0.0).is_ok());
        assert!(check_if_inbounds(&xa, 1.0).is_ok());
        assert!(check_if_inbounds(&xa, 2.0).is_ok());
        assert!(matches!(
            check_if_inbounds(&xa, 3.0).unwrap_err(),
            Domain1dError { x: 3.0 }
        ));
    }

    #[test]
    fn test_diff() {
        let s = [0.0, 1.0, -2.0, 3.0];
        let sdiff = diff(&s);

        let expected = [1.0, -3.0, 5.0];
        for i in 0..sdiff.len() {
            assert_relative_eq!(sdiff[i], expected[i]);
        }
    }
}
