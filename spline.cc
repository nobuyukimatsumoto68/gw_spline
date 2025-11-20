#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


// int main(int argc, char* argv[]){
int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------

  const int len=20;
  double x[len], y[len];

  // set data
  for (int i = 0; i < len; i++){
    x[i] = i + 0.5 * sin (i);
    y[i] = i + cos (i * i);
    printf ("%g %g\n", x[i], y[i]);
  }

  {
    // interpolate
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_spline_init (spline, x, y, len);

    // evaluation
    for (double xi = x[0]; xi < x[len]; xi += 0.01){
      double yi = gsl_spline_eval (spline, xi, acc);
      // yi = gsl_spline_eval_deriv (spline, xi, acc);
      // yi = gsl_spline_eval_deriv2 (spline, xi, acc);
      printf ("%g %g\n", xi, yi);
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  return 0;
}
