#include "grammage.h"

Grammage::Grammage() {
}

Grammage::Grammage(const PID& pid_, const Params& par) {
	pid = pid_;
	mu_c_o_2_v_A = par.mu.get() * c_light / 2. / par.v_A.get();
	H = par.galaxy_H.get();
	v_A = par.v_A.get();
	v_A_H = v_A * H;
	R_b = par.R_b.get();
	delta_hi = par.delta_hi.get();
	delta_low = par.delta_low.get();
	s = par.delta_s.get();
	D_0 = par.D_0.get();
	rho_0 = par.mu.get() / (2 * par.h_gas.get());
}

Grammage::~Grammage() {
}

double Grammage::D(const double& beta_, const double& E) const {
	double pc = pid.get_A() * beta_ * (E + proton_mass_c2);
	double R = pc / pid.get_Z();
	double x = R / R_b;
	double value = D_0 * beta_ * std::pow(x, delta_low);
	value *= std::pow(0.5 * (1. + std::pow(x, 1. / s)), (delta_hi - delta_low) * s);
	return value;
}

double Grammage::diffusion_escape_time(const double& E) const {
	return pow2(H) / D(beta(E), E);
}

double Grammage::advection_escape_time() const {
	return H / v_A;
}

double Grammage::get(const double& E) const {
	double beta_ = beta(E);
	double value = mu_c_o_2_v_A * beta_;
	value *= 1. - std::exp(-v_A_H / D(beta_, E));
	return value;
}

