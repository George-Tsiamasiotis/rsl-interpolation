#[derive(thiserror::Error, Debug)]
/// The error type for Interpolator creation and data checking.
pub enum InterpolationError {
    /// x points dataset is not sorted.
    #[error("x values must be strictly increasing.")]
    UnsortedDataset,

    /// x and y datasets have differnet length.
    #[error("Supplied datasets must be 1D and of equal length.")]
    DatasetMismatch,

    /// Supplied array size is less than the interpolation type's minimum size.
    #[error("Supplied array size is less than the interpolation type's minimum size.")]
    NotEnoughPoints,

    /// BLAS error solving Tridiagonal linear system.
    #[error("Blass error solving Tridiagonal matrix of {which_interp} Interpolator: {source}")]
    BLASTridiagError {
        which_interp: String,
        #[source]
        source: ndarray_linalg::error::LinalgError,
    },
}

#[derive(thiserror::Error, Debug)]
#[error("Supplied value is outside the range of the supplied xdata.")]
#[non_exhaustive]
/// Returned  when the supplied value is outside the range of the supplied xdata.
pub struct DomainError;
