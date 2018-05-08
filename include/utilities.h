#ifndef INCLUDE_UTILITIES_H_
#define INCLUDE_UTILITIES_H_

#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>

#include "mks.h"

double beta(const double& E);

double Gamma_Integrand(double x, void * params);

double Gamma_Integral(double slope);

#endif /* INCLUDE_UTILITIES_H_ */
