//! For 2d interpolation profiling.

use std::hint::black_box;
use std::time::Instant;

use ndarray::Array1;
use rsl_interpolation::*;

const N: usize = 100_000_000;

// ===============================================================================================

pub fn main() -> Result<(), InterpolationError> {
    let acc = &mut Accelerator::new();
    let xa: Vec<f64> = (0..100).map(f64::from).collect();
    let ya: Vec<f64> = (0..100).map(f64::from).collect();
    let interp = CubicInterpolator::build(&xa, &ya)?;

    let xs = Array1::linspace(1.0, 99.0, N);

    let start = Instant::now();
    for i in 0..N {
        let x = black_box(xs[i]);
        black_box(interp.eval(&xa, &ya, x, acc)?);
        // acc.reset();
        black_box(interp.eval_deriv(&xa, &ya, x, acc)?);
        // acc.reset();
        black_box(interp.eval_deriv(&xa, &ya, x, acc)?);
        // acc.reset();
    }
    let end = start.elapsed();
    println!("{end:?}");
    Ok(())
}
