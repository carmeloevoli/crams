#include "inelastic.h"
#include "utilities.h"

InelasticXsec::InelasticXsec() {
}

InelasticXsec::InelasticXsec(const PID& pid) {
	A = pid.get_A();
	Z = pid.get_Z();
}

InelasticXsec::~InelasticXsec() {
	std::cout << "delete inelastic sigma for particle " << A << " " << Z << "\n";
}

double InelasticXsec::get(const double& T) const {
	double sigma = (Z == 1) ? sigma_pp(T) : sigma_ST(T);
	return std::max(sigma, 1e-10 * cgs::mbarn);
}

double InelasticXsec::sigma_pp(const double& T) const {
	double value = 0;
	double x = T / E_threshold;
	if (x > 1) {
		value = 30.7 - 0.96 * log(x) + 0.18 * pow2(log(x));
		value *= pow3(1 - pow(x, -1.9));
	}
	return value * cgs::mbarn;
}

double InelasticXsec::sigma_ST(const double& T) const {
	double value = 45. * std::pow((double) A, 0.7);
	value *= 1. + 0.016 * std::sin(5.3 - 2.63 * std::log(A));
	double T_MeV = T / cgs::MeV;
	value *= 1. - 0.62 * std::exp(-T_MeV / 200.) * std::sin(10.9 * pow(T_MeV, -0.28));
	return value * cgs::mbarn;
}
