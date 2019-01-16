#ifndef INCLUDE_UTILITIES_H_
#define INCLUDE_UTILITIES_H_

#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>

#include "cgs.h"

#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define pow4(A) ((A)*(A)*(A)*(A))

double gamma_func(const double& T);
double beta_func(const double& T);
double pc_func(const int& A, const double& T);
double Gamma_Integral(double slope);

#endif /* INCLUDE_UTILITIES_H_ */
