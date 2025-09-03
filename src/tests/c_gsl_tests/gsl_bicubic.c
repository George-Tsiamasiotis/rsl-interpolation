#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline2d.h>

int main() {
  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  const size_t N = 100; /* number of points to interpolate */
  const double xa[] = {
      0.,  0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
      0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575,
      0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875,
      0.9, 0.925, 0.95, 0.975, 1.};
  const double ya[] = {
      0.,  0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
      0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575,
      0.6, 0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875,
      0.9, 0.925, 0.95, 0.975, 1.};
  const size_t nx = sizeof(xa) / sizeof(double); /* x grid points */
  const size_t ny = sizeof(ya) / sizeof(double); /* y grid points */
  double *za = malloc(nx * ny * sizeof(double));
  gsl_spline2d *spline = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  size_t i = 0;
  size_t j = 0;
  double x, y;

  /* set z grid values */

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      x = xa[i];
      y = ya[j];
      double z = cos(x) * sin(x * y);
      gsl_spline2d_set(spline, za, i, j, z);
    }
  }

  /* initialize interpolation */
  gsl_spline2d_init(spline, xa, ya, za, nx, ny);

  /* interpolate N values in x and y and print out grid for plotting */
  for (i = 0; i < N; ++i) {
    double xi = i / (N - 1.0);

    for (j = 0; j < N; ++j) {
      double yj = j / (N - 1.0);
      double z = gsl_spline2d_eval(spline, xi, yj, xacc, yacc);
      double zx = gsl_spline2d_eval_deriv_x(spline, xi, yj, xacc, yacc);
      double zy = gsl_spline2d_eval_deriv_x(spline, xi, yj, xacc, yacc);
      double zxx = gsl_spline2d_eval_deriv_xx(spline, xi, yj, xacc, yacc);
      double zyy = gsl_spline2d_eval_deriv_yy(spline, xi, yj, xacc, yacc);
      double zxy = gsl_spline2d_eval_deriv_xy(spline, xi, yj, xacc, yacc);

      printf("%f %f %f %f %f %f %f %f\n", xi, yj, z, zx, zy, zxx, zyy, zxy);
    }
    printf("\n");
  }

  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  free(za);

  return 0;
}
