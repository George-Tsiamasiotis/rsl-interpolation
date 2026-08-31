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

    /// Supplied z-grid dataset must be 1D with length of `xsize*ysize`.
    #[error("Supplied z-grid dataset must be 1D with length of `xsize*ysize`.")]
    ZGridMismatch,

    /// BLAS error solving Tridiagonal linear system.
    #[error("BLAS error solving Tridiagonal matrix of {which_interp} Interpolator: {source}")]
    BLASTridiagError {
        which_interp: String,
        #[source]
        source: ndarray_linalg::error::LinalgError,
    },

    /// Supplied value is outside the range of the supplied xdata.
    #[error("{0}")]
    Domain1dError(#[from] Domain1dError),

    /// Supplied value is outside the range of the supplied xdata or ydata.
    #[error("{0}")]
    Domain2dError(#[from] Domain2dError),
}

/// Returned when the supplied value is outside the range of the supplied xdata.
#[derive(thiserror::Error, Debug)]
#[error("1D Interpolation domain error: x = {x}")]
pub struct Domain1dError {
    pub x: f64,
}

/// Returned  when the supplied value is outside the range of the supplied xdata or ydata.
#[derive(thiserror::Error, Debug)]
#[error("1D Interpolation domain error: x = {x}, y = {y}")]
pub struct Domain2dError {
    pub x: f64,
    pub y: f64,
}
