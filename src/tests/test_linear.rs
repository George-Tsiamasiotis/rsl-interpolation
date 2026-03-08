use crate::tests::XYTable;
use crate::tests::test_interp;
use crate::{Accelerator, InterpType, Interpolation, Linear};

#[test]
fn test_type_fields() {
    let _ = <Linear as InterpType<f64>>::name(&Linear);
    let _ = <Linear as InterpType<f64>>::min_size(&Linear);
}

#[test]
fn gsl_test_linear() {
    let xa = [0.0, 1.0, 2.0, 3.0];
    let ya = [0.0, 1.0, 2.0, 3.0];

    let interp = Linear.build(&xa, &ya).unwrap();

    let xtest = [0.0, 0.5, 1.0, 1.5, 2.5, 3.0];
    let ytest = [0.0, 0.5, 1.0, 1.5, 2.5, 3.0];
    let dytest = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let iytest = [0.0, 0.125, 0.5, 9.0 / 8.0, 25.0 / 8.0, 9.0 / 2.0];

    let data_table = XYTable { x: &xa, y: &ya };

    let test_e_table = XYTable {
        x: &xtest,
        y: &ytest,
    };

    let test_d_table = XYTable {
        x: &xtest,
        y: &dytest,
    };

    let test_i_table = XYTable {
        x: &xtest,
        y: &iytest,
    };

    assert_eq!(
        interp
            .eval_deriv2(&xa, &ya, 1.5, &mut Accelerator::new())
            .unwrap(),
        0.0
    );

    test_interp(data_table, test_e_table, test_d_table, test_i_table, interp);
}
