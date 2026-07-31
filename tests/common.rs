#![allow(dead_code)]

use is_close::Comparator;
use std::f64;

use rsl_interpolation::*;

// GSL uses this to compare floats
const EPS: f64 = 1e-10;
const ATOL: f64 = 1e-9;

/// Custom Comparator with pre-set tolerances, to be used across all tests, instead of invoking the
/// macro every time.
pub(crate) fn build_comparator<'a>() -> Comparator<'a, f64> {
    let mut builder = is_close::default();
    builder.rel_tol(EPS).abs_tol(ATOL).method(is_close::AVERAGE);

    builder.compile()
}

/// A Primitive 2D table for holding the x and y values. Don't bother with num::Float here
pub(crate) struct XYTable<'a> {
    pub x: &'a [f64],
    pub y: &'a [f64],
}

/// Test function for eval(), eval_deriv() and eval_integ() for 1d interpolation. Corresponds to
/// the transferred GSL tests.
#[rustfmt::skip]
pub(crate) fn test_interp<I>(
    data_table: XYTable,
    test_e_table: XYTable,
    test_d_table: XYTable,
    test_i_table: XYTable,
    interp: I,
) where
    I: Interpolation+Interpolator,
{
    interp.name();
    interp.min_size();

    let comp = build_comparator();
    let mut acc = Accelerator::new();

    for (i, x) in test_e_table.x.iter().enumerate() {
        let s1 = interp.eval(data_table.x, data_table.y, *x, &mut acc).unwrap();
        let s2 = interp.eval_deriv(data_table.x, data_table.y, *x, &mut acc).unwrap();
        let s3 = interp.eval_integ(data_table.x, data_table.y, test_e_table.x[0], *x, &mut acc).unwrap();

        // No deriv2 tests apparently
        assert!(comp.is_close(s1, test_e_table.y[i]));
        assert!(comp.is_close(s2, test_d_table.y[i]));
        assert!(comp.is_close(s3, test_i_table.y[i]));
    }
}

/// Test function for extra tests with GSL data. Includes eval_deriv2() testing.
#[rustfmt::skip]
pub(crate) fn test_interp_extra<I>(
    data_table: XYTable,
    test_e_table: XYTable,
    test_d_table: XYTable,
    test_d2_table: XYTable,
    test_i_table: XYTable,
    interp: I,
) where
    I: Interpolation,
{
    let comp = build_comparator();
    let mut acc = Accelerator::new();

    for (i, x) in test_e_table.x.iter().enumerate() {
        let s1 = interp.eval(data_table.x, data_table.y, *x, &mut acc).unwrap();
        let s2 = interp.eval_deriv(data_table.x, data_table.y, *x, &mut acc).unwrap();
        let s3 = interp.eval_deriv2(data_table.x, data_table.y, *x, &mut acc).unwrap();
        let s4 = interp.eval_integ(data_table.x, data_table.y, test_e_table.x[0], *x, &mut acc).unwrap();

        // We need to specify an absolute tolerance, since is_close!() with abs_tol = 0 always
        // fails on 0.0, as described in https://docs.python.org/3/library/math.html#math.isclose.
        // is_close uses the python implementation.
        assert!(comp.is_close(s1, test_e_table.y[i]));
        assert!(comp.is_close(s2, test_d_table.y[i]));
        assert!(comp.is_close(s3, test_d2_table.y[i]));
        assert!(comp.is_close(s4, test_i_table.y[i]));
    }
}

// ================================================================================================

/// A Primitive 2D table for holding the x and y values. Don't bother with num::Float here
pub(crate) struct XYZTable<'a> {
    pub x: &'a [f64],
    pub y: &'a [f64],
    pub z: &'a [f64],
}

/// Test function for eval(), for 2d interpolation. Corresponds to the transferred GSL tests.
pub(crate) fn test_interp2d<I>(data_table: XYZTable, test_e_table: XYZTable, interp: I)
where
    I: Interpolation2d,
{
    let comp = build_comparator();
    let acc = &mut Accelerator2d::new();

    for (i, x) in test_e_table.x.iter().enumerate() {
        let y = test_e_table.y[i];
        let s1 = interp
            .eval(data_table.x, data_table.y, data_table.z, *x, y, acc)
            .unwrap();

        // No deriv tests apparently
        let expected = test_e_table.z[i];
        assert!(comp.is_close(s1, expected));
    }
}

/// Test function including all derivatives and iteration over all (x, y) pairs,  for use with extra
/// 2d testing.
#[rustfmt::skip]
pub(crate) fn test_interp2d_extra<I>(
    data_table: XYZTable,
    test_e_table: XYZTable,
    test_dx_table: XYZTable,
    test_dy_table: XYZTable,
    test_dxx_table: XYZTable,
    test_dyy_table: XYZTable,
    test_dxy_table: XYZTable,
    interp: I,
    interp_type: &str,
) where
    I: Interpolation2d,
{
    let comp = build_comparator();
    let acc = &mut Accelerator2d::new();

    // Access the z values linearly instead of using idx(), to comply with gsl's test output
    let mut index = 0;
    for x in test_e_table.x.iter() {
        for y in test_e_table.y.iter() {
            let eval =          interp.eval( data_table.x, data_table.y, data_table.z, *x, *y, acc).unwrap();
            let dx   =  interp.eval_deriv_x( data_table.x, data_table.y, data_table.z, *x, *y, acc).unwrap();
            let dy   =  interp.eval_deriv_y( data_table.x, data_table.y, data_table.z, *x, *y, acc).unwrap();
            let dxx  = interp.eval_deriv_xx( data_table.x, data_table.y, data_table.z, *x, *y, acc).unwrap();
            let dyy  = interp.eval_deriv_yy( data_table.x, data_table.y, data_table.z, *x, *y, acc).unwrap();
            let dxy  = interp.eval_deriv_xy( data_table.x, data_table.y, data_table.z, *x, *y, acc).unwrap();

            let expected_eval = test_e_table.z[index];
            let expected_dx= test_dx_table.z[index];
            let expected_dy= test_dy_table.z[index];
            let expected_dxx = test_dxx_table.z[index];
            let expected_dyy = test_dyy_table.z[index];
            let expected_dxy = test_dxy_table.z[index];
            index +=1;

            assert!(comp.is_close(eval, expected_eval));
            assert!(comp.is_close(dx,expected_dx));
            assert!(comp.is_close(dy,expected_dy));
            assert!(comp.is_close(dxx,expected_dxx));
            assert!(comp.is_close(dyy,expected_dyy));
            assert!(comp.is_close(dxy,expected_dxy));

        }
    }

    match interp_type {
        "bilinear" => {
            let total_evaluations = index * 4; // `deriv_xx` and `deriv_yy` do not use the cache
            assert_eq!(acc.xacc().hits() + acc.xacc().misses(), total_evaluations);
            assert_eq!(acc.yacc().hits() + acc.yacc().misses(), total_evaluations);
        }
        "bicubic" => {
            let total_evaluations = index * 6;
            assert_eq!(acc.xacc().hits() + acc.xacc().misses(), total_evaluations);
            assert_eq!(acc.yacc().hits() + acc.yacc().misses(), total_evaluations);
        }
        _ => unreachable!()
    }
}
