#include "losses.h"

Losses::Losses() {
}

Losses::Losses(const PID& pid_, const Params& par) :
		pid(pid_) {
	v_A = par.v_A.get();
	mu = par.mu.get();
	rho_0 = mu / (2 * par.h_gas.get());
}

Losses::~Losses() {
}

double Losses::adiabatic(const double& E) const {
	double value = 2 * v_A / 3 / mu / c_light;
	value *= std::sqrt(E * (E + proton_mass_c2));
	return -value;
}

double Losses::dE_dt_adiabatic(const double& E) const {
	double v = beta(E) * c_light;
	return v * adiabatic(E) * rho_0;
}

double Losses::dE_dt_ionization(const double& E) const {
	double v = beta(E) * c_light;
	return v * ionization(E) * rho_0;
}

double Losses::ionization(const double& E) const {
	double beta_ = beta(E);
	double gamma_ = (E + proton_mass_c2) / proton_mass_c2;
	double q_max = 2. * electron_mass_c2 * (pow2(gamma_) - 1.);
	q_max /= 1. + 2. * gamma_ * electron_mass_c2 / proton_mass_c2 / pid.get_A();
	double b_h = std::log(2. * electron_mass_c2 * (pow2(gamma_) - 1.) * q_max / pow2(E_H)) - 2. * pow2(beta_);
	double b_he = std::log(2. * electron_mass_c2 * (pow2(gamma_) - 1.) * q_max / pow2(E_He)) - 2. * pow2(beta_);
	double n_0 = rho_0 / proton_mass;
	double value = (2. * M_PI * pow2(electron_radius) * electron_mass_c2 * n_0) * pow2(pid.get_Z()) / beta_ * (b_h + b_he * He_abundance);
	return -value / pid.get_A() / rho_0;
}

double Losses::dE_dx(const double& E) const {
	return adiabatic(E) + ionization(E);
}

double gslLossesClassWrapper(double x, void * pp) {
	Losses * losses = (Losses *) pp;
	return losses->adiabatic(x); // + losses->ionization(x); // TODO why it does not work ionization in derivative?
}

double Losses::get_derivative(const double& E) {
	double result, abserr;

	gsl_function F;
	F.function = &gslLossesClassWrapper;
	F.params = this;
	gsl_deriv_central (&F, E, 1e-3, &result, &abserr);

	return result;
}
