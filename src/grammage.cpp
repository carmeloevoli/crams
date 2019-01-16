#include "grammage.h"

Grammage::Grammage() {
}

Grammage::Grammage(const PID& pid, const Params& params) {
	A = pid.get_A();
	Z = pid.get_Z();
	factor = params.mu * cgs::c_light / 2. / params.v_A;
	v_A = params.v_A;
	H = params.H;
	D_0 = params.D_0;
	R_b = params.R_b;
	delta = params.delta;
	ddelta = params.ddelta;
	s = params.smoothness;
}

Grammage::~Grammage() {
	std::cout << "delete grammage for particle " << A << " " << Z << "\n";
}

double Grammage::D(const double& T) const {
	double pc = pc_func(A, T);
	double R = pc / Z;
	double x = R / R_b;
	double value = beta_func(T) * std::pow(R / cgs::GeV, delta);
	value /= std::pow(1. + std::pow(x, ddelta / s), s);
	return D_0 * value;
}

double Grammage::get(const double& T) const {
	double beta = beta_func(T);
	double value = beta * factor;
	value *= 1. - std::exp(-v_A * H / D(T));
	return value;
}
