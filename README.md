## rsl-interpolation

A re-write of [`GSL's Interpolation`] in Rust.

The crates documentation can be found [`here`]

[`GSL's Interpolation`]: https://www.gnu.org/software/gsl/doc/html/interp.html
[`here`]: https://docs.rs/rsl-interpolation/latest/rsl_interpolation/

## Current status:

See [todo](TODO.md) list.

## Notes

> # **Important**
>
> In 2d Interpolation, the `za` array must be defined in **column-major (Fortran)** style. This is 
> done to comply with GSL's interface.

## Testing

All of GSL's tests have been transferred in this crate.

Additionally, some extra tests have been added, with data computed directly from GSL, to cover untested 
cases. These are located in `src/tests/*.c`. Their output can be saved and graphed with the GNU plotutils `graph` program.

GSL must be installed.


Example
```bash
make -C src/tests/c_gsl_tests akima          # Compile and run akima.c. Data is stored in src/tests/out/akima.dat
make -C src/tests/c_gsl_tests akima_plot     # Create akima.png graph
open src/tests/c_gsl_tests/plots/akima.png   # Open the graph image
make -C src/tests/c_gsl_tests clean          # Cleanup
```
