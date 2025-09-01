use crate::Accelerator;

fn setup_acc() -> Accelerator {
    Accelerator::new()
}

fn setup_xarray() -> [f64; 5] {
    [0.0, 1.0, 2.0, 3.0, 4.0]
}

#[test]
fn test_bsearch_interior_point() {
    let (acc, xarray) = (setup_acc(), setup_xarray());

    let res = acc.bsearch(&xarray, 1.5, 0, 4);
    assert_eq!(res, 1);
}

#[test]
fn test_bsearch_last_value() {
    let (acc, xarray) = (setup_acc(), setup_xarray());

    let res = acc.bsearch(&xarray, 4.0, 0, 4);
    assert_eq!(res, 3);
}

#[test]
fn test_bsearch_first_value() {
    let (acc, xarray) = (setup_acc(), setup_xarray());

    let res = acc.bsearch(&xarray, 0.0, 0, 4);
    assert_eq!(res, 0);
}

#[test]
fn test_bsearch_boundary() {
    let (acc, xarray) = (setup_acc(), setup_xarray());

    let res = acc.bsearch(&xarray, 2.0, 0, 4);
    assert_eq!(res, 2);
}

#[test]
fn test_bsearch_above_bounds() {
    let (acc, xarray) = (setup_acc(), setup_xarray());

    let res = acc.bsearch(&xarray, 10.0, 0, 4);
    assert_eq!(res, 3);
}

#[test]
fn test_bsearch_below_bounds() {
    let (acc, xarray) = (setup_acc(), setup_xarray());

    let res = acc.bsearch(&xarray, -10.0, 0, 4);
    assert_eq!(res, 0);
}

// TODO:
// #[test]
// fn bsearch_vs_lookup_bench() { }

#[test]
fn test_accelerator() {
    let xarray = setup_xarray();
    let mut acc = setup_acc();
    let mut k1 = 0;
    let mut k2 = 0;
    let mut t = false;
    let r = [
        -0.2, 0.0, 0.1, 0.7, 1.0, 1.3, 1.9, 2.0, 2.2, 2.7, 3.0, 3.1, 3.6, 4.0, 4.1, 4.9,
    ];

    // run through all the pairs of points.
    while (k1 < 16) & (k2 < 16) {
        let x = if t { r[k1] } else { r[k2] };
        t = !t;

        if !t {
            k1 = (k1 + 1) % 16;
            if k1 == 0 {
                k2 += 1;
            }
        }

        let i = acc.find(&xarray, x);
        let j = acc.bsearch(&xarray, x, 0, 4);
        assert_eq!(i, j);
    }
}
