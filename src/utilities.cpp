#include <cmath>
#include <gsl/gsl_integration.h>
#include "cgs.h"
#include "utilities.h"

#define LIMIT 1000

double beta_func(const double& T) {
	double value = std::sqrt(T * (T + cgs::proton_mass_c2));
	value /= T + cgs::proton_mass_c2;
	return value;
}

double gamma_func(const double& T) {
	double value = T + cgs::proton_mass_c2;
	value /= cgs::proton_mass_c2;
	return value;
}


double pc_func(const int& A, const double& T) {
	double value = std::sqrt(T * (T + cgs::proton_mass_c2));
	value *= (double) A;
	return value;
}

double gsl_Gamma_Integrand(double x, void * params) {
	double alpha = *(double *) params;
	double f = std::pow(x, 2. - alpha);
	f *= std::sqrt(x * x + 1) - 1;
	return f;
}

double Gamma_Integral(double slope) {
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(LIMIT);
	double result, error;
	gsl_function F;
	F.function = &gsl_Gamma_Integrand;
	F.params = &slope;
	gsl_integration_qagiu(&F, 0, 0, 1e-6, LIMIT, w, &result, &error);
	gsl_integration_workspace_free(w);
	return 4. * M_PI * result;
}

#undef LIMIT
