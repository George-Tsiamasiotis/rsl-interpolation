use crate::{DomainError, InterpolationError};

/// Checks that supplied datasets are valid.
pub(crate) fn check1d_data<T>(xa: &[T], ya: &[T], min_size: usize) -> Result<(), InterpolationError>
where
    T: crate::Num,
{
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
pub(crate) fn check2d_data<T>(
    xa: &[T],
    ya: &[T],
    za: &[T],
    min_size: usize,
) -> Result<(), InterpolationError>
where
    T: crate::Num,
{
    if (!xa.iter().is_sorted()) | (!ya.iter().is_sorted()) {
        return Err(InterpolationError::UnsortedDataset);
    }
    if (xa.len() < min_size) | (ya.len() < min_size) {
        return Err(InterpolationError::NotEnoughPoints);
    }
    if xa.len() * ya.len() != za.len() {
        Err(InterpolationError::ZGridMismatch)
    } else {
        Ok(())
    }
}

pub(crate) fn check_if_inbounds<T>(xa: &[T], x: T) -> Result<(), DomainError>
where
    T: crate::Num,
{
    if (x < *xa.first().unwrap()) | (x > *xa.last().unwrap()) {
        return Err(DomainError);
    }
    Ok(())
}

/// Calculates the n-th discrete difference: out[i] = s[i+1]-s[i].
pub(crate) fn diff<T>(s: &[T]) -> Vec<T>
where
    T: crate::Num,
{
    s.windows(2)
        .map(|xy| {
            let [x, y] = xy else { unreachable!() };
            *y - *x
        })
        .collect::<Vec<T>>()
}

/// Function for doing the spline integral evaluation, which is common to both the cspline and
/// akima methods.
pub(crate) fn integ_eval<T>(ai: T, bi: T, ci: T, di: T, xi: T, a: T, b: T) -> T
where
    T: crate::Num,
{
    let quarter = T::from(0.25).unwrap();
    let half = T::from(0.5).unwrap();
    let third = T::from(1.0 / 3.0).unwrap();

    let r1 = a - xi;
    let r2 = b - xi;
    let r12 = r1.powi(2);
    let r22 = r2.powi(2);
    let rsum = r1 + r2;
    let bterm = half * bi * rsum;
    let cterm = third * ci * (r12 + r22 + r1 * r2);
    let dterm = quarter * di * rsum * (r12 + r22);

    (b - a) * (ai + bterm + cterm + dterm)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_check_data() {
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
    fn test_check_if_inbounds() {
        let xa = [0.0, 1.0, 2.0];

        assert!(check_if_inbounds(&xa, 0.0).is_ok());
        assert!(check_if_inbounds(&xa, 1.0).is_ok());
        assert!(check_if_inbounds(&xa, 2.0).is_ok());
        assert!(matches!(
            check_if_inbounds(&xa, 3.0).unwrap_err(),
            DomainError
        ));
    }

    #[test]
    fn test_diff() {
        let s = [0.0, 1.0, -2.0, 3.0];
        let sdiff = diff(&s);

        assert_eq!(sdiff, vec![1.0, -3.0, 5.0]);
    }
}
