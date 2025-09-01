use is_close::is_close;
use std::f64;

// GSL uses this to compare floats
const EPS: f64 = 1e-10;

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

/// Test function for all Interpolation Types
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
