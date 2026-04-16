/*
 * Test for bicubic interpolator, including all the derivatives and iteration
 * over all (x, y) pairs.
 * */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline2d.h>

int main() {
  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  const size_t N = 9; /* number of points to interpolate */
  double xa[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  double ya[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
  /* least common multiple of x and y */
  double za[] = {1, 2,  3,  4,  5,  6,  7,  8,  2, 2, 6,  4,  10, 6,  14, 8,
                 3, 6,  3,  12, 15, 6,  21, 24, 4, 4, 12, 4,  20, 12, 28, 8,
                 5, 10, 15, 20, 5,  30, 35, 40, 6, 6, 6,  12, 30, 6,  42, 24,
                 7, 14, 21, 28, 35, 42, 7,  56, 8, 8, 24, 8,  40, 24, 56, 8};

  const size_t nx = sizeof(xa) / sizeof(double); /* x grid points */
  const size_t ny = sizeof(ya) / sizeof(double); /* y grid points */
  gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  size_t i = 0;
  size_t j = 0;

  /* initialize interpolation */
  gsl_spline2d_init(spline, xa, ya, za, nx, ny);

  /* interpolate N values in x and y and print out grid for plotting */
  for (i = 0; i < N; ++i) {
    double xi = xa[0] + i * (xa[nx - 1] - xa[0]) / (N - 1);

    for (j = 0; j < N; ++j) {
      double yj = ya[0] + j * (ya[ny - 1] - ya[0]) / (N - 1);
      double z = gsl_spline2d_eval(spline, xi, yj, xacc, yacc);
      double dx = gsl_spline2d_eval_deriv_x(spline, xi, yj, xacc, yacc);
      double dy = gsl_spline2d_eval_deriv_y(spline, xi, yj, xacc, yacc);
      double dxx = gsl_spline2d_eval_deriv_xx(spline, xi, yj, xacc, yacc);
      double dyy = gsl_spline2d_eval_deriv_yy(spline, xi, yj, xacc, yacc);
      double dxy = gsl_spline2d_eval_deriv_xy(spline, xi, yj, xacc, yacc);

      printf("%.15f          %.15f          %.15f          %.15f          "
             "%.15f          %.15f          %.15f         %.15f\n",
             xi, yj, z, dx, dy, dxx, dyy, dxy);
    }
    printf("\n");
  }

  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);

  return 0;
}
