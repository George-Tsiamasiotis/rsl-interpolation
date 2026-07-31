use rsl_interpolation::*;

#[test]
fn dyn_interp() {
    let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    let mut acc = Accelerator::new();

    let x = 0.5;
    let dyn_interp = DynInterpolator::build(Interpolation1dType::Cubic, &xa, &ya).unwrap();
    let a = dyn_interp.eval(&xa, &ya, x, &mut acc).unwrap();
    let _ = format!("{:?}", &dyn_interp);

    let clone = dyn_interp.clone();
    drop(dyn_interp);
    let b = clone.eval(&xa, &ya, x, &mut acc).unwrap();
    assert_eq!(a, b);
}

#[test]
fn spline() {
    let xa = [0.0, 1.0, 2.0, 3.0, 4.0];
    let ya = [0.0, 2.0, 4.0, 6.0, 8.0];
    let mut acc = Accelerator::new();

    let x = 0.5;
    let interpolator = DynInterpolator::build(Interpolation1dType::Cubic, &xa, &ya).unwrap();
    let spline = Spline::new(interpolator, &xa, &ya);
    let a = spline.eval(x, &mut acc).unwrap();
    let _ = format!("{:?}", &spline);

    let clone = spline.clone();
    drop(spline);
    let b = clone.eval(x, &mut acc).unwrap();
    assert_eq!(a, b);
}

#[test]
fn dyn_interp2d() {
    let xa = [0.0, 1.0, 2.0];
    let ya = [0.0, 2.0, 4.0];
    let za = [0.0, 1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 5.0, 6.0];
    let mut acc = Accelerator2d::new();

    let x = 0.5;
    let y = 1.0;
    let dyn_interp =
        DynInterpolator2d::build(Interpolation2dType::Bilinear, &xa, &ya, &za).unwrap();
    let a = dyn_interp.eval(&xa, &ya, &za, x, y, &mut acc).unwrap();
    let _ = format!("{:?}", &dyn_interp);

    let clone = dyn_interp.clone();
    drop(dyn_interp);
    let b = clone.eval(&xa, &ya, &za, x, y, &mut acc).unwrap();
    assert_eq!(a, b);
}

#[test]
fn spline2d() {
    let xa = [0.0, 1.0, 2.0];
    let ya = [0.0, 2.0, 4.0];
    let za = [0.0, 1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 5.0, 6.0];
    let mut acc = Accelerator2d::new();

    let x = 0.5;
    let y = 1.0;
    let interpolator =
        DynInterpolator2d::build(Interpolation2dType::Bilinear, &xa, &ya, &za).unwrap();
    let spline = Spline2d::new(interpolator, &xa, &ya, &za);
    let a = spline.eval(x, y, &mut acc).unwrap();
    let _ = format!("{:?}", &spline);

    let clone = spline.clone();
    drop(spline);
    let b = clone.eval(x, y, &mut acc).unwrap();
    assert_eq!(a, b);
}
