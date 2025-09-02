use is_close::is_close;
use std::f64;

// GSL uses this to compare floats
const EPS: f64 = 1e-10;
const ATOL: f64 = 1e-9;

use crate::Accelerator;
use crate::Interpolation;
use crate::Interpolation2d;

mod test_accel;
mod test_akima;
mod test_cubic;
mod test_cubic_periodic;
mod test_linear;
mod test_steffen;

mod test_bilinear;

/// A Primitive 2D table for holding the x and y values. Don't bother with num::Float here
pub(crate) struct XYTable<'a> {
    x: &'a [f64],
    y: &'a [f64],
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
    I: Interpolation<f64>,
{
    let mut acc = Accelerator::new();

    for (i, x) in test_e_table.x.iter().enumerate() {
        let s1 = interp.eval(data_table.x, data_table.y, *x, &mut acc).unwrap();
        let s2 = interp.eval_deriv(data_table.x, data_table.y, *x, &mut acc).unwrap();
        let s3 = interp.eval_integ(data_table.x, data_table.y, test_e_table.x[0], *x, &mut acc).unwrap();

        // No deriv2 tests apparently
        assert!(is_close!(s1, test_e_table.y[i], rel_tol = EPS));
        assert!(is_close!(s2, test_d_table.y[i], rel_tol = EPS));
        assert!(is_close!(s3, test_i_table.y[i], rel_tol = EPS));
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
    I: Interpolation<f64>,
{
    let mut acc = Accelerator::new();

    for (i, x) in test_e_table.x.iter().enumerate() {
        let s1 = interp.eval(data_table.x, data_table.y, *x, &mut acc).unwrap();
        let s2 = interp.eval_deriv(data_table.x, data_table.y, *x, &mut acc).unwrap();
        let s3 = interp.eval_deriv2(data_table.x, data_table.y, *x, &mut acc).unwrap();
        let s4 = interp.eval_integ(data_table.x, data_table.y, test_e_table.x[0], *x, &mut acc).unwrap();

        // We need to specify an absolute tolerance, since is_close!() with abs_tol = 0 always
        // fails on 0.0, as described in https://docs.python.org/3/library/math.html#math.isclose.
        // is_close uses the python implementation.
        assert!(is_close!(s1, test_e_table.y[i], rel_tol = EPS, abs_tol = ATOL));
        assert!(is_close!(s2, test_d_table.y[i], rel_tol = EPS, abs_tol = ATOL));
        assert!(is_close!(s3, test_d2_table.y[i], rel_tol = EPS, abs_tol = ATOL));
        assert!(is_close!(s4, test_i_table.y[i], rel_tol = EPS, abs_tol = ATOL));
    }
}

// ================================================================================================

/// A Primitive 2D table for holding the x and y values. Don't bother with num::Float here
pub(crate) struct XYZTable<'a> {
    x: &'a [f64],
    y: &'a [f64],
    z: &'a [f64],
}

/// Test function for eval(), for 2d interpolation. Corresponds to the transferred GSL tests.
#[rustfmt::skip]
pub(crate) fn test_interp2d<I>(data_table: XYZTable, test_e_table: XYZTable, interp: I)
where
    I: Interpolation2d<f64>,
{
    let mut xacc = Accelerator::new();
    let mut yacc = Accelerator::new();

    for (i, x) in test_e_table.x.iter().enumerate() {
        let y = test_e_table.y[i];
        let s1 = interp.eval(data_table.x, data_table.y, data_table.z, *x, y, &mut xacc, &mut yacc).unwrap();

        // No deriv tests apparently
        let expected = test_e_table.z[i];
        assert!(is_close!(s1, expected, rel_tol = EPS));
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
) where
    I: Interpolation2d<f64>,
{
    let mut xacc = Accelerator::new();
    let mut yacc = Accelerator::new();

    // Access the z values linearly instead of using idx(), to comply with gsl's test output
    let mut index = 0; 
    for x in test_e_table.x.iter() {
        for y in test_e_table.y.iter() {
            let eval =          interp.eval( data_table.x, data_table.y, data_table.z, *x, *y, &mut xacc, &mut yacc).unwrap();
            let dx   =  interp.eval_deriv_x( data_table.x, data_table.y, data_table.z, *x, *y, &mut xacc, &mut yacc).unwrap();
            let dy   =  interp.eval_deriv_y( data_table.x, data_table.y, data_table.z, *x, *y, &mut xacc, &mut yacc).unwrap();
            let dxx  = interp.eval_deriv_xx( data_table.x, data_table.y, data_table.z, *x, *y, &mut xacc, &mut yacc).unwrap();
            let dyy  = interp.eval_deriv_yy( data_table.x, data_table.y, data_table.z, *x, *y, &mut xacc, &mut yacc).unwrap();
            let dxy  = interp.eval_deriv_xy( data_table.x, data_table.y, data_table.z, *x, *y, &mut xacc, &mut yacc).unwrap();

            let expected_eval = test_e_table.z[index];
            let expected_eval_deriv_x = test_dx_table.z[index];
            let expected_eval_deriv_y = test_dy_table.z[index];
            let expected_eval_deriv_xx = test_dxx_table.z[index];
            let expected_eval_deriv_yy = test_dyy_table.z[index];
            let expected_eval_deriv_xy = test_dxy_table.z[index];
            index +=1;

            assert!(is_close!(eval, expected_eval, rel_tol = EPS, abs_tol = ATOL));
            assert!(is_close!(dx, expected_eval_deriv_x, rel_tol = EPS, abs_tol = ATOL));
            assert!(is_close!(dy, expected_eval_deriv_y, rel_tol = EPS, abs_tol = ATOL));
            assert!(is_close!(dxx, expected_eval_deriv_xx, rel_tol = EPS));
            assert!(is_close!(dyy, expected_eval_deriv_yy, rel_tol = EPS));
            assert!(is_close!(dxy, expected_eval_deriv_xy, rel_tol = EPS));
        }
    }
}
