## rsl-interpolation

A re-write of [`GSL's Interpolation`] in Rust.

The crates documentation can be found [`here`]

[`GSL's Interpolation`]: https://www.gnu.org/software/gsl/doc/html/interp.html
[`here`]: https://docs.rs/rsl-interpolation/latest/rsl_interpolation/

## Current status:

See [todo](TODO.md) list.

## Testing

All of GSL's tests have been transferred in this crate.

Additionally, some extra tests have been added, with data computed directly from GSL, to cover untested 
cases. These are located in `src/gsl_tests/*.c`. Their output can be saved and graphed with the GNU plotutils `graph` program.

GSL must be installed.


Example
```bash
make -C src/tests akima        # Compile and run akima.c. Data is stored in src/tests/out/akima.dat
make -C src/tests akima_plot   # Create akima.png graph
open src/tests/out/akima.png   # Open the graph image
make -C src/tests clean        # Cleanup
```
