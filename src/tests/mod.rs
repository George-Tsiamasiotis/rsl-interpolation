use is_close::is_close;
use std::f64;

// GSL uses this to compare floats
const EPS: f64 = 1e-10;
const ATOL: f64 = 1e-9;

use crate::Accelerator;
use crate::Interpolation;

mod test_accel;
mod test_akima;
mod test_cubic;
mod test_cubic_periodic;
mod test_linear;
mod test_steffen;

/// A Primitive 2D table for holding the x and y values. Don't bother with num::Float here
pub(crate) struct XYTable<'a> {
    x: &'a [f64],
    y: &'a [f64],
}

/// Test function for eval(), eval_deriv() and eval_integ(). Corresponds to the transferred GSL
/// tests.
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
        let s1 = interp
            .eval(data_table.x, data_table.y, *x, &mut acc)
            .unwrap();
        let s2 = interp
            .eval_deriv(data_table.x, data_table.y, *x, &mut acc)
            .unwrap();
        let s3 = interp
            .eval_integ(data_table.x, data_table.y, test_e_table.x[0], *x, &mut acc)
            .unwrap();

        // No deriv2 tests apparently
        assert!(is_close!(s1, test_e_table.y[i], rel_tol = EPS));
        assert!(is_close!(s2, test_d_table.y[i], rel_tol = EPS));
        assert!(is_close!(s3, test_i_table.y[i], rel_tol = EPS));
    }
}

/// Test function for extra tests with GSL data. Includes eval_deriv2() testing.
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
        let s1 = interp
            .eval(data_table.x, data_table.y, *x, &mut acc)
            .unwrap();
        let s2 = interp
            .eval_deriv(data_table.x, data_table.y, *x, &mut acc)
            .unwrap();
        let s3 = interp
            .eval_deriv2(data_table.x, data_table.y, *x, &mut acc)
            .unwrap();
        let s4 = interp
            .eval_integ(data_table.x, data_table.y, test_e_table.x[0], *x, &mut acc)
            .unwrap();

        // We need to specify an absolute tolerance, since is_close!() with abs_tol = 0 always
        // fails on 0.0, as described in https://docs.python.org/3/library/math.html#math.isclose.
        // is_close uses the python implementation.
        #[rustfmt::skip]
        assert!(is_close!( s1, test_e_table.y[i], rel_tol = EPS, abs_tol = ATOL));
        #[rustfmt::skip]
        assert!(is_close!( s2, test_d_table.y[i], rel_tol = EPS, abs_tol = ATOL));
        #[rustfmt::skip]
        assert!(is_close!( s3, test_d2_table.y[i], rel_tol = EPS, abs_tol = ATOL));
        #[rustfmt::skip]
        assert!(is_close!( s4, test_i_table.y[i], rel_tol = EPS, abs_tol = ATOL));
    }
}
