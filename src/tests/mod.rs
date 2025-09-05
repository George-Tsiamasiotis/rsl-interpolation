use is_close::Comparator;
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

mod test_bicubic;
mod test_bilinear;

/// Custom Comparator with pre-set tolerances, to be used across all tests, instead of invoking the
/// macro every time.
fn build_comparator<'a, T>() -> Comparator<'a, T>
where
    T: crate::Num + 'a,
{
    let mut builder = is_close::default::<T>();
    builder
        .rel_tol(T::from(EPS).unwrap())
        .abs_tol(T::from(ATOL).unwrap())
        .method(is_close::AVERAGE);

    builder.compile()
}

/// A Primitive 2D table for holding the x and y values. Don't bother with num::Float here
pub(crate) struct XYTable<'a, T>
where
    T: crate::Num,
{
    x: &'a [T],
    y: &'a [T],
}

/// Test function for eval(), eval_deriv() and eval_integ() for 1d interpolation. Corresponds to
/// the transferred GSL tests.
#[rustfmt::skip]
pub(crate) fn test_interp<I, T>(
    data_table: XYTable<T>,
    test_e_table: XYTable<T>,
    test_d_table: XYTable<T>,
    test_i_table: XYTable<T>,
    interp: I,
) where
    T: crate::Num,
    I: Interpolation<T>,
{
    let comp = build_comparator::<T>();
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
pub(crate) fn test_interp_extra<I, T>(
    data_table: XYTable<T>,
    test_e_table: XYTable<T>,
    test_d_table: XYTable<T>,
    test_d2_table: XYTable<T>,
    test_i_table: XYTable<T>,
    interp: I,
) where
    T: crate::Num,
    I: Interpolation<T>,
{
    let comp = build_comparator::<T>();
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
pub(crate) struct XYZTable<'a, T> {
    x: &'a [T],
    y: &'a [T],
    z: &'a [T],
}

/// Test function for eval(), for 2d interpolation. Corresponds to the transferred GSL tests.
#[rustfmt::skip]
pub(crate) fn test_interp2d<I, T>(data_table: XYZTable<T>, test_e_table: XYZTable<T>, interp: I)
where
    T: crate::Num,
    I: Interpolation2d<T>,
{
    let comp = build_comparator::<T>();
    let mut xacc = Accelerator::new();
    let mut yacc = Accelerator::new();

    for (i, x) in test_e_table.x.iter().enumerate() {
        let y = test_e_table.y[i];
        let s1 = interp.eval(data_table.x, data_table.y, data_table.z, *x, y, &mut xacc, &mut yacc).unwrap();

        // No deriv tests apparently
        let expected = test_e_table.z[i];
        assert!(comp.is_close(s1, expected));
    }
}

/// Test function including all derivatives and iteration over all (x, y) pairs,  for use with extra 
/// 2d testing.
#[rustfmt::skip]
pub(crate) fn test_interp2d_extra<I, T>(
    data_table: XYZTable<T>,
    test_e_table: XYZTable<T>,
    test_dx_table: XYZTable<T>,
    test_dy_table: XYZTable<T>,
    test_dxx_table: XYZTable<T>,
    test_dyy_table: XYZTable<T>,
    test_dxy_table: XYZTable<T>,
    interp: I,
) where
    T: crate::Num,
    I: Interpolation2d<T>,
{
    let comp = build_comparator::<T>();
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
}
