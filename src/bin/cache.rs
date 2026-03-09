//! For 2d interpolation profiling

use std::hint::black_box;
use std::time::Instant;

use ndarray::Array1;
use rsl_interpolation::*;

const N: usize = 100_000_000;
const TYP: Bilinear = Bilinear; // keep to "Bilinear", otherwise needs BLAS to build and publish

// ===============================================================================================

pub fn main() {
    let xacc = &mut Accelerator::new();
    let yacc = &mut Accelerator::new();
    let cache = &mut Cache::new();
    let xa: Vec<f64> = (0..100).map(f64::from).collect();
    let ya: Vec<f64> = (0..200).map(f64::from).collect();
    let za: Vec<f64> = (0..(100 * 200)).map(f64::from).collect();
    let interp2d = TYP.build(&xa, &ya, &za).unwrap();

    let xs = Array1::linspace(1.0, 99.0, N);
    let ys = Array1::linspace(1.0, 199.0, N);

    let start = Instant::now();
    for i in 0..N {
        let x = black_box(xs[i]);
        let y = black_box(ys[i]);
        black_box(
            interp2d
                .eval(&xa, &ya, &za, x, y, xacc, yacc, cache)
                .unwrap(),
        );
        // clear_cache(xacc, yacc, cache);
        black_box(
            interp2d
                .eval_deriv_x(&xa, &ya, &za, x, y, xacc, yacc, cache)
                .unwrap(),
        );
        // clear_cache(xacc, yacc, cache);
        black_box(
            interp2d
                .eval_deriv_y(&xa, &ya, &za, x, y, xacc, yacc, cache)
                .unwrap(),
        );
        // clear_cache(xacc, yacc, cache);
    }
    let end = start.elapsed();
    println!("{end:?}")
}

#[allow(dead_code)]
fn clear_cache(xacc: &mut Accelerator, yacc: &mut Accelerator, cache: &mut Cache<f64>) {
    xacc.reset();
    yacc.reset();
    cache.soft_reset();
}
