#include <cmath>
#include <iostream>
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

std::vector<double> LinAxis(const double& min, const double& max, const size_t& size) {
	const double dx = (max - min) / (double) (size - 1);
	std::vector<double> v(size);
	for (size_t i = 0; i < size; ++i) {
		double value = min + dx * i;
		v[i] = value;
	}
	return v;
}

std::vector<double> LogAxis(const double& min, const double& max, const size_t& size) {
	const double delta_log = std::exp(std::log(max / min) / (size - 1));
	std::vector<double> v(size);
	for (size_t i = 0; i < size; ++i) {
		double value = std::exp(std::log(min) + (double) i * std::log(delta_log));
		v[i] = value;
	}
	return v;
}

double LinearInterpolator(const std::vector<double>& x, const std::vector<double>& y,
		const double& x_new) {
	auto value = double();
	if (x_new > x.front() && x_new <= x.back()) {
		size_t const i = std::lower_bound(x.begin(), x.end(), x_new) - x.begin();
		double t = (x_new - x.at(i - 1)) / (x.at(i) - x.at(i - 1));
		value = y.at(i - 1) * (1. - t) + y.at(i) * t;
	}
	return value;
}

//double LinearInterpolatorLog(const std::vector<double>& x, const std::vector<double>& y,
//		const double& x_new) {
//	auto value = double();
//	if (x_new >= x.front() && x_new <= x.back()) {
//		size_t const i = std::lower_bound(x.begin(), x.end(), x_new) - x.begin();
//		double t = std::log(x_new) - std::log(x.at(i - 1));
//		t /= std::log(x.at(i)) - std::log(x.at(i - 1));
//		value = std::log(y.at(i - 1)) * (1. - t) + std::log(y.at(i)) * t;
//		value = std::exp(value);
//	}
//	return value;
//}

size_t getLowerIndex(const std::vector<double>& v, double x) {
	assert(x >= v.front());
	size_t i = 0;
	while (v.at(i + 1) < x)
		i++;
	return i;
}

double LinearInterpolatorLog(const std::vector<double>& x, const std::vector<double>& y,
		const double& x_new) {
	if (x_new < x.front() || x_new > x.back())
		return 0;
	else {
		size_t const i = getLowerIndex(x, x_new);
		double t = std::log(x_new) - std::log(x.at(i));
		t /= std::log(x.at(i + 1)) - std::log(x.at(i));
		double v = std::log(y.at(i)) * (1. - t) + std::log(y.at(i + 1)) * t;
		return std::exp(v);
	}
}

#undef LIMIT
