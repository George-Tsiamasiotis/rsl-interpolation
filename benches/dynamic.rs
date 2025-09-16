//! Benchmark of static vs dynamic Spline

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use rsl_interpolation::*;

const WARMUP_MILLIS: u64 = 100;
const MEASUREMENT_SECS: u64 = 3;
const XA: [f64; 10] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
const YA: [f64; 10] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
const X: f64 = 3.0;
const TYP: Akima = Akima;

// ===============================================================================================

pub fn static_vs_dynamic_spline1d(c: &mut Criterion) {
    let mut acc = Accelerator::new();
    let spline = Spline::new(TYP, &XA, &YA).unwrap();
    let dyn_spline = Spline::new_dyn(TYP, &XA, &YA).unwrap();

    let mut group = c.benchmark_group("Spline vs DynSpline");
    group.warm_up_time(Duration::from_millis(WARMUP_MILLIS));
    group.measurement_time(Duration::from_secs(MEASUREMENT_SECS));

    group.bench_function("Spline", |b| b.iter(|| spline.eval(X, &mut acc)));
    group.bench_function("DynSpline", |b| b.iter(|| dyn_spline.eval(X, &mut acc)));

    group.finish();
}

criterion_group!(benches, static_vs_dynamic_spline1d);
criterion_main!(benches);
