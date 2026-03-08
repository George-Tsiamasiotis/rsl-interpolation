## rsl-interpolation

A re-write of [`GSL's Interpolation`] in Rust.

The crates documentation can be found [`here`]

[`GSL's Interpolation`]: https://www.gnu.org/software/gsl/doc/html/interp.html
[`here`]: https://docs.rs/rsl-interpolation/latest/rsl_interpolation/

## Current status:

See [todo](TODO.md) list.

## Notes

1. `rsl-interpolation` requires LAPACK FFI, so you must use **just one** of the corresponding [`ndarray_linalg features`](https://github.com/rust-ndarray/ndarray-linalg?tab=readme-ov-file#backend-features).

2. In 2d Interpolation, the `za` array must be defined in **column-major (Fortran)** style. This is done to comply with GSL's interface.

## Testing

All of GSL's tests have been transferred in this crate.

Additionally, some extra tests have been added, with data computed directly from GSL, to cover untested cases. These are located in `src/tests/c_gsl_tests/*.c`. Their output can be saved and graphed with the GNU plotutils `graph` program.

`GSL` must be installed.
For the plots to work, `gnuplot` must be installed.

Example
```bash
make -C src/tests/c_gsl_tests                # Run tests and create plots
open src/tests/c_gsl_tests/plots/akima.png   # Open the graph image
make -C src/tests/c_gsl_tests clean          # Cleanup
```
