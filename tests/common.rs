#![allow(dead_code)]

use approx::assert_relative_eq;
use rsl_interpolation::*;

// GSL uses this to compare floats
pub const EPS: f64 = 1e-10;
pub const ATOL: f64 = 1e-9;

/// A Primitive 2D table for holding the x and y values. Don't bother with num::Float here
pub(crate) struct XYTable<'a> {
    pub x: &'a [f64],
    pub y: &'a [f64],
}

/// Test function for eval(), eval_deriv() and eval_integ() for 1d interpolation. Corresponds to
/// the transferred GSL tests.
pub(crate) fn test_interp(
    test_e_table: XYTable,
    test_d_table: XYTable,
    test_i_table: XYTable,
    spline: Spline,
) {
    let mut acc = Accelerator::new();

    for (i, x) in test_e_table.x.iter().enumerate() {
        let s1 = spline.eval(*x, &mut acc).unwrap();
        let s2 = spline.eval_deriv(*x, &mut acc).unwrap();
        let s3 = spline.eval_integ(test_e_table.x[0], *x, &mut acc).unwrap();

        // No deriv2 tests apparently
        assert_relative_eq!(s1, test_e_table.y[i], epsilon = EPS);
        assert_relative_eq!(s2, test_d_table.y[i], epsilon = EPS);
        assert_relative_eq!(s3, test_i_table.y[i], epsilon = EPS);
    }
}

/// Test function for extra tests with GSL data. Includes eval_deriv2() testing.
pub(crate) fn test_interp_extra(
    test_e_table: XYTable,
    test_d_table: XYTable,
    test_d2_table: XYTable,
    test_i_table: XYTable,
    spline: Spline,
) {
    let mut acc = Accelerator::new();

    for (i, x) in test_e_table.x.iter().enumerate() {
        let s1 = spline.eval(*x, &mut acc).unwrap();
        let s2 = spline.eval_deriv(*x, &mut acc).unwrap();
        let s3 = spline.eval_deriv2(*x, &mut acc).unwrap();
        let s4 = spline.eval_integ(test_e_table.x[0], *x, &mut acc).unwrap();

        assert_relative_eq!(s1, test_e_table.y[i], epsilon = EPS);
        assert_relative_eq!(s2, test_d_table.y[i], epsilon = EPS);
        assert_relative_eq!(s3, test_d2_table.y[i], epsilon = EPS);
        assert_relative_eq!(s4, test_i_table.y[i], epsilon = EPS);
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
pub(crate) fn test_interp2d(test_e_table: XYZTable, spline: Spline2d) {
    let acc = &mut Accelerator2d::new();

    for (i, x) in test_e_table.x.iter().enumerate() {
        let y = test_e_table.y[i];
        let s1 = spline.eval(*x, y, acc).unwrap();

        // No deriv tests apparently
        let expected = test_e_table.z[i];
        assert_relative_eq!(s1, expected, epsilon = EPS);
    }
}

/// Test function including all derivatives and iteration over all (x, y) pairs,  for use with extra
/// 2d testing.
pub(crate) fn test_interp2d_extra(
    test_e_table: XYZTable,
    test_dx_table: XYZTable,
    test_dy_table: XYZTable,
    test_dxx_table: XYZTable,
    test_dyy_table: XYZTable,
    test_dxy_table: XYZTable,
    spline: Spline2d,
    interp_type: &str,
) {
    let acc = &mut Accelerator2d::new();

    // Access the z values linearly instead of using idx(), to comply with gsl's test output
    let mut index = 0;
    for x in test_e_table.x.iter() {
        for y in test_e_table.y.iter() {
            let eval = spline.eval(*x, *y, acc).unwrap();
            let dx = spline.eval_deriv_x(*x, *y, acc).unwrap();
            let dy = spline.eval_deriv_y(*x, *y, acc).unwrap();
            let dxx = spline.eval_deriv_xx(*x, *y, acc).unwrap();
            let dyy = spline.eval_deriv_yy(*x, *y, acc).unwrap();
            let dxy = spline.eval_deriv_xy(*x, *y, acc).unwrap();

            let expected_eval = test_e_table.z[index];
            let expected_dx = test_dx_table.z[index];
            let expected_dy = test_dy_table.z[index];
            let expected_dxx = test_dxx_table.z[index];
            let expected_dyy = test_dyy_table.z[index];
            let expected_dxy = test_dxy_table.z[index];
            index += 1;

            assert_relative_eq!(eval, expected_eval, epsilon = EPS);
            assert_relative_eq!(dx, expected_dx, epsilon = EPS);
            assert_relative_eq!(dy, expected_dy, epsilon = EPS);
            assert_relative_eq!(dxx, expected_dxx, epsilon = EPS);
            assert_relative_eq!(dyy, expected_dyy, epsilon = EPS);
            assert_relative_eq!(dxy, expected_dxy, epsilon = EPS);
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
        _ => unreachable!(),
    }
}
