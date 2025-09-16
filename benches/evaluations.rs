//! General benchmark of low-level interpolation evaluation methods.

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use rsl_interpolation::*;

const WARMUP_MILLIS: u64 = 100;
const MEASUREMENT_SECS: u64 = 3;
const XA: [f64; 10] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
const YA: [f64; 10] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
const X: f64 = 3.0;
const A: f64 = 1.5;
const B: f64 = 6.8;
const TYP: Cubic = Cubic;

// ===============================================================================================

pub fn interp_evals(c: &mut Criterion) {
    let mut acc = Accelerator::new();
    let interp = TYP.build(&XA, &YA).unwrap();

    let mut group = c.benchmark_group("Low-level Interpolator Evaluation Methods");
    group.warm_up_time(Duration::from_millis(WARMUP_MILLIS));
    group.measurement_time(Duration::from_secs(MEASUREMENT_SECS));

    group.bench_function("interp.eval()", |b| {
        b.iter(|| interp.eval(&XA, &YA, X, &mut acc))
    });
    group.bench_function("interp.eval_deriv()", |b| {
        b.iter(|| interp.eval_deriv(&XA, &YA, X, &mut acc))
    });
    group.bench_function("interp.eval_deriv2()", |b| {
        b.iter(|| interp.eval_deriv2(&XA, &YA, X, &mut acc))
    });
    group.bench_function("interp.eval_integ()", |b| {
        b.iter(|| interp.eval_integ(&XA, &YA, A, B, &mut acc))
    });
}

criterion_group!(benches, interp_evals);
criterion_main!(benches);
