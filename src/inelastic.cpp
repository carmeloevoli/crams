#include "inelastic.h"

InelasticXsec::InelasticXsec() {
}

InelasticXsec::InelasticXsec(const PID& pid_) :
		pid(pid_) {
}

InelasticXsec::~InelasticXsec() {
}

double InelasticXsec::get(const double& E) const {
	double sigma = (pid.is_H()) ? sigma_pp(E) : sigma_ST(E, pid.get_A());
	return sigma;
}

double InelasticXsec::sigma_pp(const double& E) const {
	double value = 0;
	double x = E / E_threshold;
	if (x > 1) {
		value = 30.7 - 0.96 * log(x) + 0.18 * pow2(log(x));
		value *= pow3(1 - pow(x, -1.9));
	}
	return value * mbarn;
}

double InelasticXsec::sigma_ST(const double& E, const int& A_k) const {
	double value = 45. * std::pow((double) A_k, 0.7);
	value *= 1. + 0.016 * std::sin(5.3 - 2.63 * std::log(A_k));
	value *= 1. - 0.62 * std::exp(-E / (0.2 * GeV)) * std::sin(10.9 * pow(E / MeV, -0.28));
	return value * mbarn;
}

