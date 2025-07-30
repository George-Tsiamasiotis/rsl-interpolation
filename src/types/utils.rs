use crate::{DomainError, InterpolationError};

/// Checks that supplied datasets are valid.
pub(crate) fn check_data<T>(xa: &[T], ya: &[T], min_size: usize) -> Result<(), InterpolationError>
where
    T: num::Float,
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

pub(crate) fn check_if_inbounds<T>(xa: &[T], x: T) -> Result<(), DomainError>
where
    T: num::Float,
{
    if (x < *xa.first().unwrap()) | (x > *xa.last().unwrap()) {
        return Err(DomainError);
    }
    Ok(())
}

/// Calculates the n-th discrete difference: out[i] = s[i+1]-s[i].
pub(crate) fn diff<T>(s: &[T]) -> Vec<T>
where
    T: num::Float,
{
    s.windows(2)
        .map(|xy| {
            let [x, y] = xy else { unreachable!() };
            *y - *x
        })
        .collect::<Vec<T>>()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_check_data() {
        let ya = [0.0, 1.0, 2.0];

        let xa = [0.0, 1.0, 2.0];
        assert!(check_data(&xa, &ya, 2).is_ok());
        assert!(matches!(
            check_data(&xa, &ya, 4).unwrap_err(),
            InterpolationError::NotEnoughPoints
        ));

        let xa = [2.0, 1.0, 2.0];
        assert!(matches!(
            check_data(&xa, &ya, 2).unwrap_err(),
            InterpolationError::UnsortedDataset
        ));

        let xa = [0.0, 1.0, 2.0, 3.0];
        assert!(matches!(
            check_data(&xa, &ya, 2).unwrap_err(),
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
