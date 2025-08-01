use crate::Akima;
use crate::Interpolation;
use crate::gsl_tests::XYTable;
use crate::gsl_tests::test_interp;

#[test]
fn test_akima() {
    let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    let ya = [0.0, 1.0, 2.0, 3.0, 4.0];

    let xtest = [0.0, 0.5, 1.0, 2.0];
    let ytest = [0.0, 0.5, 1.0, 2.0];
    let dytest = [1.0, 1.0, 1.0, 1.0];
    let iytest = [0.0, 0.125, 0.5, 2.0];

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

    let interp = Akima::new(&xa, &ya).unwrap();
    test_interp(data_table, test_e_table, test_d_table, test_i_table, interp);
}
