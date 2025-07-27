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
