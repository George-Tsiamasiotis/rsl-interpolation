//! General benchmark of low-level 2d interpolation evaluation methods.

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use rsl_interpolation::*;

const WARMUP_MILLIS: u64 = 100;
const MEASUREMENT_SECS: u64 = 3;
const XA: [f64; 10] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
const YA: [f64; 10] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
const X: f64 = 3.4;
const Y: f64 = 7.1;
const TYP: Bicubic = Bicubic;

// ===============================================================================================

pub fn interp2d_evals(c: &mut Criterion) {
    let mut xacc = Accelerator::new();
    let mut yacc = Accelerator::new();
    let za: Vec<f64> = (0..100).map(f64::from).collect();
    let interp2d = TYP.build(&XA, &YA, &za).unwrap();

    let mut group = c.benchmark_group("Low-level 2D Interpolator Evaluation Methods");
    group.warm_up_time(Duration::from_millis(WARMUP_MILLIS));
    group.measurement_time(Duration::from_secs(MEASUREMENT_SECS));

    group.bench_function("interp2d.eval()", |b| {
        b.iter(|| interp2d.eval(&XA, &YA, &za, X, Y, &mut xacc, &mut yacc))
    });
    group.bench_function("interp2d.eval_deriv_x()", |b| {
        b.iter(|| interp2d.eval_deriv_x(&XA, &YA, &za, X, Y, &mut xacc, &mut yacc))
    });
    group.bench_function("interp2d.eval_deriv_y()", |b| {
        b.iter(|| interp2d.eval_deriv_y(&XA, &YA, &za, X, Y, &mut xacc, &mut yacc))
    });
    group.bench_function("interp2d.eval_deriv_xx()", |b| {
        b.iter(|| interp2d.eval_deriv_xx(&XA, &YA, &za, X, Y, &mut xacc, &mut yacc))
    });
    group.bench_function("interp2d.eval_deriv_yy()", |b| {
        b.iter(|| interp2d.eval_deriv_yy(&XA, &YA, &za, X, Y, &mut xacc, &mut yacc))
    });
    group.bench_function("interp2d.eval_deriv_xy()", |b| {
        b.iter(|| interp2d.eval_deriv_xy(&XA, &YA, &za, X, Y, &mut xacc, &mut yacc))
    });
}

criterion_group!(benches, interp2d_evals);
criterion_main!(benches);
