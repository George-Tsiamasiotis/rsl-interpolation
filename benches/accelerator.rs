//! Accelerator benchmark
//!
//! Use Bicubic, since it is by far the most computationally expensive
//! evaulation, and resetting the Accelerator should be negligible.

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ndarray::Array1;
use rsl_interpolation::*;

const WARMUP_MILLIS: u64 = 100;
const MEASUREMENT_SECS: u64 = 3;

// ===============================================================================================

pub fn interp_accel_eval(c: &mut Criterion) {
    let mut acc = Accelerator::new();
    let xa: Vec<f64> = (0..301).map(f64::from).collect();
    let ya: Vec<f64> = (0..301).map(f64::from).collect();
    let interp = Linear.build(&xa, &ya).unwrap();

    let mut group = c.benchmark_group("1D Interpolation Accelerator (Incremental values)");
    group.warm_up_time(Duration::from_millis(WARMUP_MILLIS));
    group.measurement_time(Duration::from_secs(MEASUREMENT_SECS));

    let arr = Array1::linspace(1.0, 300.0, 10000).to_vec();

    group.bench_function("1D eval with Accelerator", |b| {
        b.iter(|| {
            for x in arr.iter() {
                interp.eval(&xa, &ya, *x, &mut acc).unwrap();
            }
        })
    });
    // println!("{acc:#?}");

    group.bench_function("1D eval without Accelerator", |b| {
        b.iter(|| {
            for x in arr.iter() {
                acc.reset();
                interp.eval(&xa, &ya, *x, &mut acc).unwrap();
            }
        })
    });

    group.finish();
}

pub fn interp2d_accel_eval(c: &mut Criterion) {
    let mut xacc = Accelerator::new();
    let mut yacc = Accelerator::new();
    let mut cache = Cache::new();
    let xa: Vec<f64> = (0..101).map(f64::from).collect();
    let ya: Vec<f64> = (0..101).map(f64::from).collect();
    let za: Vec<f64> = (0..10201).map(f64::from).collect();
    let interp2d = Bicubic.build(&xa, &ya, &za).unwrap();

    let mut group = c.benchmark_group("2D Interpolation Accelerator");
    group.warm_up_time(Duration::from_millis(WARMUP_MILLIS));
    group.measurement_time(Duration::from_secs(MEASUREMENT_SECS));

    let xarr = Array1::linspace(1.0, 10.0, 100).to_vec();
    let yarr = Array1::linspace(1.0, 10.0, 100).to_vec();

    group.bench_function("2D eval with Accelerator", |b| {
        b.iter(|| {
            for x in xarr.iter() {
                for y in yarr.iter() {
                    interp2d
                        .eval(&xa, &ya, &za, *x, *y, &mut xacc, &mut yacc, &mut cache)
                        .unwrap();
                }
            }
        })
    });

    // println!("{xacc:#?}");
    // println!("{yacc:#?}");

    group.bench_function("2D eval without Accelerator", |b| {
        b.iter(|| {
            for x in xarr.iter() {
                xacc.reset();
                for y in yarr.iter() {
                    yacc.reset();
                    interp2d
                        .eval(&xa, &ya, &za, *x, *y, &mut xacc, &mut yacc, &mut cache)
                        .unwrap();
                }
            }
        })
    });

    group.finish();
}

criterion_group!(benches, interp_accel_eval, interp2d_accel_eval);
criterion_main!(benches);
