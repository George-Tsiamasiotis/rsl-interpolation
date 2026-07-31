//! Tests for the correct Cache-Accelerator hits/misses synchronization fix.

mod common;

use common::*;
use rsl_interpolation::*;

#[test]
#[rustfmt::skip]
fn test_cache_accelerator_update() {
    let xa = [0.0, 1.0, 2.0, 3.0];
    let ya = [0.0, 1.0, 2.0, 3.0];
    #[rustfmt::skip]
    let za = [
        1.0, 1.1, 1.2, 1.3,
        1.1, 1.2, 1.3, 1.4,
        1.2, 1.3, 1.4, 1.5,
        1.3, 1.4, 1.5, 1.6,
    ];

    let interp = BicubicInterpolator::build(&xa, &ya, &za).unwrap();
    let acc = &mut Accelerator2d::new();

    interp.eval(&xa, &ya, &za, 0.6, 0.6, acc).unwrap(); // hit, hit

    assert_eq!(acc.xacc().hits(), 1);
    assert_eq!(acc.yacc().hits(), 1);
    assert_eq!(acc.xacc().misses(), 0);
    assert_eq!(acc.yacc().misses(), 0);

    let acc = &mut Accelerator2d::new();

    interp.eval(&xa, &ya, &za, 2.5, 2.5, acc).unwrap(); // miss, miss
    interp.eval(&xa, &ya, &za, 0.5, 0.5, acc).unwrap(); // miss, miss
    interp.eval(&xa, &ya, &za, 0.6, 0.6, acc).unwrap(); // hit, hit
    interp.eval(&xa, &ya, &za, 0.7, 0.7, acc).unwrap(); // hit, hit
    interp.eval(&xa, &ya, &za, 1.5, 1.5, acc).unwrap(); // miss, miss
    interp.eval(&xa, &ya, &za, 1.6, 1.6, acc).unwrap(); // hit, hit
    interp.eval(&xa, &ya, &za, 2.5, 2.5, acc).unwrap(); // miss, miss
    interp.eval(&xa, &ya, &za, 0.5, 0.5, acc).unwrap(); // miss, miss

    assert_eq!(acc.xacc().hits(), 3);
    assert_eq!(acc.yacc().hits(), 3);
    assert_eq!(acc.xacc().misses(), 5);
    assert_eq!(acc.yacc().misses(), 5);

}

/// Tests for an erroneous cache update when calling `deriv_x` and `deriv_y` with an uninitialized
/// Cache, that came up with the Cache-Accelerator synchronization fix. Bilinear first derivatives
/// are the only ones using the cache.
#[test]
#[rustfmt::skip]
fn test_cache_accelerator_update_bilinear_first_deriv() {
    // Results from GSL's tests.
    let xa = [0.0, 1.0, 2.0, 3.0];
    let ya = [0.0, 1.0, 2.0, 3.0];
    #[rustfmt::skip]
    let za = [
        1.0, 1.1, 1.2, 1.4,
        1.3, 1.4, 1.5, 1.7,
        1.5, 1.6, 1.7, 1.9,
        1.6, 1.9, 2.2, 2.3,
    ];

    let interp = BilinearInterpolator::build(&xa, &ya, &za).unwrap();
    let acc = &mut Accelerator2d::new();

    let comp = build_comparator();

    let dzdy = interp.eval(&xa, &ya, &za, 0.5, 0.5, acc).unwrap();
    assert!(comp.is_close(dzdy, 1.2));
    assert_eq!(acc.xacc().hits(), 1);
    assert_eq!(acc.yacc().hits(), 1);
    assert_eq!(acc.xacc().misses() , 0);
    assert_eq!(acc.yacc().misses() , 0);
    assert_eq!(acc.xacc().cache(), 0);
    assert_eq!(acc.yacc().cache(), 0);

    let dzdy = interp.eval_deriv_y(&xa, &ya, &za, 1.5, 1.5, acc).unwrap();
    assert!(comp.is_close(dzdy, 0.2));
    assert_eq!(acc.xacc().hits(), 1);
    assert_eq!(acc.yacc().hits(), 1);
    assert_eq!(acc.xacc().misses() , 1);
    assert_eq!(acc.yacc().misses() , 1);
    assert_eq!(acc.xacc().cache(), 1);
    assert_eq!(acc.yacc().cache(), 1);

    let dzdy = interp.eval_deriv_y(&xa, &ya, &za, 1.5, 3.0, acc).unwrap();
    assert!(comp.is_close(dzdy, 0.4));
    assert_eq!(acc.xacc().hits(), 2);
    assert_eq!(acc.yacc().hits(), 1);
    assert_eq!(acc.xacc().misses() , 1);
    assert_eq!(acc.yacc().misses() , 2);
    assert_eq!(acc.xacc().cache(), 1);
    assert_eq!(acc.yacc().cache(), 2);
}
