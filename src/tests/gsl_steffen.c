/* Test for Steffen interpolation, with f(x) = 1 + x^2
 * */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

int main(void) {

  int i;
  double x[6] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  size_t data_points = 6;
  size_t test_points = 19;
  double xi, yi, y[data_points];

  printf("#m=0,S=17\n");

  for (i = 0; i < data_points; i++) {
    y[i] = 1 + pow(x[i], 2);
    printf("%.15f\t %.15f\n", x[i], y[i]);
  }

  printf("#m=1,S=0\n");

  {
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, data_points);

    gsl_spline_init(spline, x, y, data_points);

    double step = x[4] / (double)test_points;

    for (xi = x[0]; xi <= x[data_points - 1]; xi += step) {
      yi = gsl_spline_eval(spline, xi, acc);
      printf("%.15f\t %.15f\n", xi, yi);
    }

    printf("#m=2,S=0\n");

    for (xi = x[0]; xi <= x[data_points - 1]; xi += step) {
      yi = gsl_spline_eval_deriv(spline, xi, acc);
      printf("%.15f\t %.15f\n", xi, yi);
    }

    printf("#m=3,S=0\n");

    for (xi = x[0]; xi <= x[data_points - 1]; xi += step) {
      yi = gsl_spline_eval_deriv2(spline, xi, acc);
      printf("%.15f\t %.15f\n", xi, yi);
    }

    printf("#m=4,S=0\n");

    for (xi = x[0]; xi <= x[data_points - 1]; xi += step) {
      yi = gsl_spline_eval_integ(spline, x[0], xi, acc);
      printf("%.15f\t %.15f\n", xi, yi);
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  return 0;
}
