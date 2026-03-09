//! Tests for the correct Cache-Accelerator hits/misses synchronization fix.

use crate::tests::*;
use crate::*;

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

    let interp = Bicubic.build(&xa, &ya, &za).unwrap();

    let xacc = &mut Accelerator::new();
    let yacc = &mut Accelerator::new();
    let cache = &mut Cache::new();

    interp.eval(&xa, &ya, &za, 0.6, 0.6, xacc, yacc, cache).unwrap(); // hit, hit

    assert_eq!(xacc.hits, 1);
    assert_eq!(yacc.hits, 1);
    assert_eq!(xacc.misses, 0);
    assert_eq!(yacc.misses, 0);

    let xacc = &mut Accelerator::new();
    let yacc = &mut Accelerator::new();
    let cache = &mut Cache::new();

    interp.eval(&xa, &ya, &za, 2.5, 2.5, xacc, yacc, cache).unwrap(); // miss, miss
    interp.eval(&xa, &ya, &za, 0.5, 0.5, xacc, yacc, cache).unwrap(); // miss, miss
    interp.eval(&xa, &ya, &za, 0.6, 0.6, xacc, yacc, cache).unwrap(); // hit, hit
    interp.eval(&xa, &ya, &za, 0.7, 0.7, xacc, yacc, cache).unwrap(); // hit, hit
    interp.eval(&xa, &ya, &za, 1.5, 1.5, xacc, yacc, cache).unwrap(); // miss, miss
    interp.eval(&xa, &ya, &za, 1.6, 1.6, xacc, yacc, cache).unwrap(); // hit, hit
    interp.eval(&xa, &ya, &za, 2.5, 2.5, xacc, yacc, cache).unwrap(); // miss, miss
    interp.eval(&xa, &ya, &za, 0.5, 0.5, xacc, yacc, cache).unwrap(); // miss, miss

    assert_eq!(xacc.hits, 3);
    assert_eq!(yacc.hits, 3);
    assert_eq!(xacc.misses, 5);
    assert_eq!(yacc.misses, 5);

}

/// Tests for an erroneous cache update when calling `deriv_x` and `deriv_y` with an uninitialized
/// Cache, that came up with the Cache-Accelerator synchronization fix. Bilinear first derivatives
/// are the only
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

    let interp = Bilinear.build(&xa, &ya, &za).unwrap();
    let xacc = &mut Accelerator::new();
    let yacc = &mut Accelerator::new();
    let cache = &mut Cache::new();

    let comp = build_comparator();

    let dzdy = interp.eval(&xa, &ya, &za, 0.5, 0.5, xacc, yacc, cache).unwrap();
    assert!(comp.is_close(dzdy, 1.2));
    assert_eq!(xacc.hits, 1);
    assert_eq!(yacc.hits, 1);
    assert_eq!(xacc.misses , 0);
    assert_eq!(yacc.misses , 0);
    assert_eq!(xacc.cache, 0);
    assert_eq!(yacc.cache, 0);
    assert_eq!(cache.get_xy_indices(), (0, 0));

    let dzdy = interp.eval_deriv_y(&xa, &ya, &za, 1.5, 1.5, xacc, yacc, cache).unwrap();
    assert!(comp.is_close(dzdy, 0.2));
    assert_eq!(xacc.hits, 1);
    assert_eq!(yacc.hits, 1);
    assert_eq!(xacc.misses , 1);
    assert_eq!(yacc.misses , 1);
    assert_eq!(xacc.cache, 1);
    assert_eq!(yacc.cache, 1);
    assert_eq!(cache.get_xy_indices(), (1, 1));

    let dzdy = interp.eval_deriv_y(&xa, &ya, &za, 1.5, 3.0, xacc, yacc, cache).unwrap();
    assert!(comp.is_close(dzdy, 0.4));
    assert_eq!(xacc.hits, 2);
    assert_eq!(yacc.hits, 1);
    assert_eq!(xacc.misses , 1);
    assert_eq!(yacc.misses , 2);
    assert_eq!(xacc.cache, 1);
    assert_eq!(yacc.cache, 2);
    assert_eq!(cache.get_xy_indices(), (1, 2));
}
