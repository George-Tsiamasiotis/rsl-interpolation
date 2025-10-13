use std::time::Instant;

use rsl_interpolation::*;

const X: f64 = 3.4;
const Y: f64 = 7.1;
const TYP: Bilinear = Bilinear;

// ===============================================================================================

pub fn main() {
    let mut xacc = Accelerator::new();
    let mut yacc = Accelerator::new();
    let mut cache = Cache::new();
    let xa: Vec<f64> = (0..100).map(f64::from).collect();
    let ya: Vec<f64> = (0..3620).map(f64::from).collect();
    let za: Vec<f64> = (0..(100 * 3620)).map(f64::from).collect();
    let interp2d = TYP.build(&xa, &ya, &za).unwrap();

    let start = Instant::now();
    for _ in 0..10000000 {
        interp2d
            .eval(&xa, &ya, &za, X, Y, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        interp2d
            .eval_deriv_x(&xa, &ya, &za, X, Y, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
        interp2d
            .eval_deriv_y(&xa, &ya, &za, X, Y, &mut xacc, &mut yacc, &mut cache)
            .unwrap();
    }
    let end = start.elapsed();
    println!("{end:?}")
}
