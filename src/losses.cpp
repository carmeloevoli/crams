#include <cmath>
#include <string>
#include <gsl/gsl_deriv.h>
#include "cgs.h"
#include "losses.h"
#include "utilities.h"

Losses::Losses() {
}

Losses::Losses(const PID& pid, const Params& params) {
	A = pid.get_A();
	Z = pid.get_Z();
	factor_ad = 2. * params.v_A / 3. / params.mu / cgs::c_light;
	mu = params.mu;
}

Losses::~Losses() {
	std::cout << "delete losses for particle " << A << " " << Z << "\n";

}

double Losses::get(const double& T) const {
	return dE_dx_adiabatic(T) + dE_dx_ionization(T);
}

double Losses::dE_dx_adiabatic(const double& T) const {
	double value = factor_ad * std::sqrt(T * (T + 2. * cgs::proton_mass_c2));
	return -value;
}

double Losses::dE_dx_ionization(const double& T) const {
	double beta = beta_func(T);
	double gamma = gamma_func(T);
	constexpr double me_c2 = cgs::electron_mass_c2;
	constexpr double two_PI_re2 = 2. * M_PI * pow2(cgs::electron_radius);
	double Q_max = 2. * cgs::electron_mass_c2 * pow2(beta) * pow2(gamma);
	Q_max /= 1. + 2. * gamma * cgs::electron_mass / ((double) A * cgs::proton_mass);
	double B_H = std::log(2. * me_c2 * (pow2(gamma) - 1.) * Q_max / pow2(cgs::Is_H));
	B_H -= 2. * pow2(beta);
	double B_He = std::log(2. * me_c2 * (pow2(gamma) - 1.) * Q_max / pow2(cgs::Is_He));
	B_He -= 2. * pow2(beta);
	double dEdx = two_PI_re2 * me_c2 * (double) pow2(Z) * (B_H + B_He * cgs::f_He);
	dEdx /= cgs::proton_mass * (1. + 4. * cgs::f_He) * (double) A * pow2(beta);
	return -dEdx;
}

/* double Losses::dE_dt_adiabatic(const double& T) const {
 double v = beta_func(T) * cgs::c_light;
 return v * dE_dx_adiabatic(T) * cgs::rho_ism;
 }*/

/* double Losses::dE_dt_ionization(const double& T) const {
 double v = beta_func(T) * cgs::c_light;
 return v * dE_dx_ionization(T) * cgs::rho_ism;
 }*/

double gslLossesClassWrapper(double x, void * pp) {
	Losses * losses = (Losses *) pp;
	return losses->dE_dx_adiabatic(x) + losses->dE_dx_ionization(x);
}

double Losses::get_derivative(const double& T) {
	double result, abserr;
	gsl_function F;
	F.function = &gslLossesClassWrapper;
	F.params = this;
	gsl_deriv_central(&F, T, 1e-5, &result, &abserr);
	return result;
}
