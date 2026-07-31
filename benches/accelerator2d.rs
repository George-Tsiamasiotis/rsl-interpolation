//! Accelerator benchmark

use std::time::Duration;

use criterion::{Criterion, criterion_group, criterion_main};
use ndarray::{Array1, Array2};
use rand::{RngExt, SeedableRng};
use rsl_interpolation::*;

const WARMUP_MILLIS: u64 = 100;
const MEASUREMENT_SECS: u64 = 3;
const EVAL_ARRAY_LEN: usize = 10_000;

// ===============================================================================================

pub fn interp_accel_eval(c: &mut Criterion) {
    let mut acc = Accelerator2d::new();
    let xa: Vec<f64> = (0..301).map(f64::from).collect();
    let ya: Vec<f64> = (0..301).map(f64::from).collect();
    let za: Array2<f64> = Array2::from_shape_fn((301, 301), |(i, j)| xa[i] * ya[j]);
    let za = za.flatten_with_order(ndarray::Order::ColumnMajor).to_vec();
    let interp = BicubicInterpolator::build(&xa, &ya, &za).unwrap();

    let mut group = c.benchmark_group("2D Interpolation Accelerator");
    group.warm_up_time(Duration::from_millis(WARMUP_MILLIS));
    group.measurement_time(Duration::from_secs(MEASUREMENT_SECS));

    // ====================================

    let arr = Array1::from_elem(EVAL_ARRAY_LEN, 4.0);

    group.bench_function("2D eval with Accelerator (Same value)", |b| {
        b.iter(|| {
            for x in arr.iter() {
                interp.eval(&xa, &ya, &za, *x, *x, &mut acc).unwrap();
            }
        })
    });
    println!(
        "Same value: {:#?}\n{:#?}",
        acc.xacc().clone(),
        acc.yacc().clone()
    );

    group.bench_function("2D eval without Accelerator (Same value)", |b| {
        b.iter(|| {
            for x in arr.iter() {
                acc.reset();
                interp.eval(&xa, &ya, &za, *x, *x, &mut acc).unwrap();
            }
        })
    });

    // ====================================

    let arr = Array1::linspace(0.0, 300.0, EVAL_ARRAY_LEN);

    group.bench_function("2D eval with Accelerator (Incremental values)", |b| {
        b.iter(|| {
            for x in arr.iter() {
                interp.eval(&xa, &ya, &za, *x, *x, &mut acc).unwrap();
            }
        })
    });
    println!(
        "Incremental values: {:#?}\n{:#?}",
        acc.xacc().clone(),
        acc.yacc().clone()
    );

    group.bench_function("2D eval without Accelerator (Incremental values)", |b| {
        b.iter(|| {
            for x in arr.iter() {
                acc.reset();
                interp.eval(&xa, &ya, &za, *x, *x, &mut acc).unwrap();
            }
        })
    });

    // ====================================

    let mut rng = rand::rngs::SmallRng::seed_from_u64(42);
    let arr: [f64; EVAL_ARRAY_LEN] = core::array::from_fn(|_| rng.random_range(0.0..300.0));

    group.bench_function("2D eval with Accelerator (Random values)", |b| {
        b.iter(|| {
            for x in arr.iter() {
                interp.eval(&xa, &ya, &za, *x, *x, &mut acc).unwrap();
            }
        })
    });
    println!(
        "Random values: {:#?}\n{:#?}",
        acc.xacc().clone(),
        acc.yacc().clone()
    );

    group.bench_function("2D eval without Accelerator (Random values)", |b| {
        b.iter(|| {
            for x in arr.iter() {
                acc.reset();
                interp.eval(&xa, &ya, &za, *x, *x, &mut acc).unwrap();
            }
        })
    });

    group.finish();
}

criterion_group!(benches, interp_accel_eval);
criterion_main!(benches);
