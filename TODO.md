## GSL features to be implemented

- [ ] 1D Interpolation
	- [x] Instantiation [`gsl_interp_alloc()`, `gsl_interp_init()`, `gsl_interp_free()`]
	- [x] 1D Interpolation types [`gsl_interp_type`]
		- [x] Linear [`gsl_interp_linear`]
		- [ ] Polynomial [`gsl_interp_polynomial`] [^1]
		- [x] Cubic [`gsl_interp_cspline`]
		- [x] Cubic Periodic [`gsl_interp_cspline_periodic`] **only works for 3 points at the moment; the general case is missing a cyclically tridiagonal matrix solver, which is currently not implemented by [`ndarray_linalg`]**.
		- [x] Akima [`gsl_interp_akima`], could use some better testing
		- [x] Akima Periodic [`gsl_interp_akima_periodic`]
		- [x] Steffen [`gsl_interp_steffen`]
	- [x] Evaluation [^2]
		- [x] f(x) evaluation [`gsl_interp_eval()`]
		- [x] f'(x) evaluation [`gsl_interp_eval_deriv()`]
 		- [x] f''(x) evaluation [`gsl_interp_eval_deriv2()`]
 		- [x] Numerical Integral [`gsl_interp_integ()`]
	- [x] Utility functions
		- [x] Name [`gsl_interp_name()`]
		- [x] Minimum number of points [`gsl_interp_min_size()` and `gsl_interp_type_min_size()`]
	- [ ] Higher level Interface (Splines)

---

- [ ] 2D Interpolation
	- [x] Instantiation [`gsl_interp2d_alloc()`, `gsl_interp2d_init()`, `gsl_interp2d_free()`]
	- [x] 2D Interpolation Grids [`gsl_interp2d_set`, `gsl_interp2d_get()`, `gsl_interp2d_idx`]
	- [ ] 1D Interpolation types [`gsl_interp2d_type`]
		- [x] Bilinear [`gsl_interp2d_bilinear`]
		- [ ] Bicubic [`gsl_interp2d_bicubic`]
	- [x] Utility functions
		- [x] Name [`gsl_interp2d_name()`]
		- [x] Minimum number of points [`gsl_interp2d_min_size()` and `gsl_interp2d_type_min_size()`]
	- [x] Evaluation [^2]
		- [x] f(x, y) evaluation [`gsl_interp2d_eval()`]
 		- [x] f(x, y) extrapolated evaluation [`gsl_interp2d_eval_extrap()`]
		- [x] fx(x, y) evaluation [`gsl_interp2d_eval_deriv_x()`]
		- [x] fy(x, y) evaluation [`gsl_interp2d_eval_deriv_y()`]
		- [x] fxx(x, y) evaluation [`gsl_interp2d_eval_deriv_xx()`]
		- [x] fyy(x, y) evaluation [`gsl_interp2d_eval_deriv_yy()`]
		- [x] fxy(x, y) evaluation [`gsl_interp2d_eval_deriv_xy()`]
	- [ ] Higher level Interface (Splines)
	
---

- [x] Acceleration
	- [x] Instantiation [`gsl_interp_accel_alloc()`, `gsl_interp_accel_reset()`, `gsl_interp_accel_free()`]
	- [x] Lookup [`gsl_interp_bsearch()`, `gsl_interp_accel_find()`]

[`ndarray_linalg`]: https://docs.rs/ndarray-linalg/latest/

[^1]: Polynomial Interpolation is probably not going to be implemented, since it relies on divided differences and after all, it is a terrible interpolation method.

[^2]: `_e()` evaluation functions are probably not gonna be implemented.
