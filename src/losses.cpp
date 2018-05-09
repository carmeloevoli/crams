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

double Losses::ionization(const double& T) const {
	double E_tot = double(pid.get_A()) * (T + proton_mass_c2);
	double M_A = proton_mass_c2 * double(pid.get_A());
	double p_c = std::sqrt(pow2(E_tot) - pow2(M_A));
	double gamma_ = (T + proton_mass_c2) / proton_mass_c2;
	double beta_ = p_c / E_tot;
	double n_0 = rho_0 / proton_mass;
	double IsH = 19. * eV;
	double IsHe = 44. * eV;
	double Q_max = 2 * electron_mass_c2 * pow2(beta_) * pow2(gamma_) / (1 + 2 * gamma_ * electron_mass_c2 / M_A); // GeV
	double B_H = std::log(2 * electron_mass_c2 * (pow2(gamma_) - 1) * Q_max / IsH / IsH) - 2 * pow2(beta_);
	double B_He = std::log(2 * electron_mass_c2 * (pow2(gamma_) - 1) * Q_max / IsHe / IsHe) - 2 * pow2(beta_);
	double dEdx_ion = (2 * M_PI * pow2(electron_radius) * electron_mass_c2 * n_0) * (double(pow2(pid.get_Z())) / beta_) * (B_H + B_He * He_abundance);
	dEdx_ion /= rho_0 / double(pid.get_A());
	return -dEdx_ion;
}

double Losses::dE_dx(const double& E) const {
	return adiabatic(E) + ionization(E);
}

double gslLossesClassWrapper(double x, void * pp) {
	Losses * losses = (Losses *) pp;
	return losses->adiabatic(x) + losses->ionization(x);
}

double Losses::get_derivative(const double& E) {
	double result, abserr;

	gsl_function F;
	F.function = &gslLossesClassWrapper;
	F.params = this;

	gsl_deriv_forward(&F, E, 1e-3, &result, &abserr);

	return result;
}
