use std::{f64::consts::PI, time::Duration};

use criterion::{Criterion, criterion_group, criterion_main};
use ndarray::Array1;
use rsl_interpolation::*;

const WARMUP_MILLIS: u64 = 100;
const MEASUREMENT_SECS: u64 = 3;

// ===============================================================================================

pub fn interp2d_cache(c: &mut Criterion) {
    let mut xacc = Accelerator::new();
    let mut yacc = Accelerator::new();
    let mut cache = Cache::new();
    let xa: Vec<f64> = (0..101).map(f64::from).collect();
    let ya: Vec<f64> = (0..3620).map(f64::from).collect();
    let za: Vec<f64> = (0..101 * 3620).map(f64::from).collect();
    let interp2d = Bicubic.build(&xa, &ya, &za).unwrap();

    let mut group = c.benchmark_group("2D Interpolation Cache");
    group.warm_up_time(Duration::from_millis(WARMUP_MILLIS));
    group.measurement_time(Duration::from_secs(MEASUREMENT_SECS));

    let xarr = Array1::linspace(0.0, 0.1, 100).to_vec();
    let yarr = Array1::linspace(0.0, PI, 100).to_vec();

    group.bench_function("2D eval with Cache", |b| {
        b.iter(|| {
            for i in 0..xarr.len() {
                let x = xarr[i];
                let y = yarr[i];
                interp2d
                    .eval(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                interp2d
                    .eval_deriv_x(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                interp2d
                    .eval_deriv_y(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                interp2d
                    .eval_deriv_xx(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                interp2d
                    .eval_deriv_yy(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                interp2d
                    .eval_deriv_xy(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
            }
        })
    });

    println!("{xacc:?}");
    println!("{yacc:?}");

    group.bench_function("2D eval without Cache", |b| {
        b.iter(|| {
            for i in 0..xarr.len() {
                let x = xarr[i];
                let y = yarr[i];
                interp2d
                    .eval(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                cache.reset();
                interp2d
                    .eval_deriv_x(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                cache.reset();
                interp2d
                    .eval_deriv_y(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                cache.reset();
                interp2d
                    .eval_deriv_xx(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                cache.reset();
                interp2d
                    .eval_deriv_yy(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                cache.reset();
                interp2d
                    .eval_deriv_xy(&xa, &ya, &za, x, y, &mut xacc, &mut yacc, &mut cache)
                    .unwrap();
                cache.reset();
            }
        })
    });

    group.finish();
}

criterion_group!(benches, interp2d_cache);
criterion_main!(benches);
