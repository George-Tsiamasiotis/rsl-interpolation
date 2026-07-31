use rsl_interpolation::*;

#[test]
fn test_z_idx() {
    // [
    //      0, 0, 0,
    //      0, 0, 1,
    //      0, 0, 0,
    //      0, 0, 0,
    // ]
    let shape = (4, 3); // Fortran style
    assert_eq!(z_idx(1, 2, shape.0, shape.1), 9);
}

#[test]
#[should_panic]
fn test_z_idx_panic() {
    let shape = (4, 3);
    z_idx(10, 200, shape.0, shape.1);
}

#[test]
fn test_z_set() {
    let xa = [0.0, 1.0];
    let ya = [0.0, 2.0];

    #[rustfmt::skip]
    let mut za = [
        0.0, 1.0,
        1.0, 0.5,
    ];

    let za00 = 100.0;
    let za01 = 300.0;
    let za10 = 200.0;
    let za11 = 400.0;

    let xlen = xa.len();
    let ylen = ya.len();

    z_set(&mut za, za00, 0, 0, xlen, ylen);
    z_set(&mut za, za01, 0, 1, xlen, ylen);
    z_set(&mut za, za10, 1, 0, xlen, ylen);
    z_set(&mut za, za11, 1, 1, xlen, ylen);

    assert_eq!(za, [100.0, 200.0, 300.0, 400.0,]);
}

#[test]
#[should_panic]
fn test_z_set_panic() {
    let xa = [0.0, 1.0];
    let ya = [0.0, 2.0];

    #[rustfmt::skip]
    let mut za = [
        0.0, 1.0,
        1.0, 0.5,
    ];

    let xlen = xa.len();
    let ylen = ya.len();

    z_set(&mut za, 1.0, 10, 10000, xlen, ylen);
}

#[test]
fn test_z_get() {
    let xa = [0.0, 1.0, 2.0];
    let ya = [0.0, 1.0, 2.0, 3.0];
    #[rustfmt::skip]
    let za = [
        0.0, 1.0, 2.0,
        3.0, 4.0, 5.0,
        6.0, 5.0, 4.0,
        3.0, 99.0, 1.0, // we want 99.0
    ];

    let (i, j) = (1, 3);
    let idx = z_get(&za, i, j, xa.len(), ya.len());
    let expected = 99.0;
    assert_eq!(idx, expected);
}

#[test]
#[should_panic]
fn test_z_get_panic() {
    let xa = [0.0, 1.0];
    let ya = [0.0, 1.0, 3.0];
    #[rustfmt::skip]
    let za = [
        0.0, 1.0,
        3.0, 4.0,
        6.0, 5.0,
    ];
    z_get(&za, 10, 2000, xa.len(), ya.len());
}
