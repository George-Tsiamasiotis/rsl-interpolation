const DOMAIN_ERROR_MSG: &str =
    "Supplied value is outside the range of the supplied xdata or ydata.";

#[derive(thiserror::Error, Debug)]
/// The error type for Interpolator creation and data checking.
pub enum InterpolationError {
    /// x points dataset is not sorted.
    #[error("x values must be strictly increasing.")]
    UnsortedDataset,

    /// x and y datasets have different length.
    #[error("Supplied datasets must be 1D and of equal length.")]
    DatasetMismatch,

    /// Supplied array size is less than the interpolation type's minimum size.
    #[error("Supplied array size is less than the interpolation type's minimum size.")]
    NotEnoughPoints,

    /// Suppled z-grid dataset must be 1D with length of `xsize*ysize`.
    #[error("Suppled z-grid dataset must be 1D with length of `xsize*ysize`.")]
    ZGridMismatch,

    /// BLAS error solving Tridiagonal linear system.
    #[error("BLAS error solving Tridiagonal matrix of {which_interp} Interpolator: {source}")]
    BLASTridiagError {
        which_interp: String,
        #[source]
        source: ndarray_linalg::error::LinalgError,
    },

    /// Supplied value is outside the range of the supplied xdata or ydata.
    #[error("{DOMAIN_ERROR_MSG}")]
    DomainError(#[from] DomainError),

    /// Invalid Interpolation Type
    #[error("`{0}`: Invalid Interpolation Type")]
    InvalidType(Box<str>),
}

#[derive(thiserror::Error, Debug)]
#[error("{DOMAIN_ERROR_MSG}")]
#[non_exhaustive]
/// Returned  when the supplied value is outside the range of the supplied xdata or ydata.
pub struct DomainError;
