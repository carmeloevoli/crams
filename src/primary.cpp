#include "primary.h"

PrimarySource::PrimarySource() {
}

PrimarySource::PrimarySource(const PID& pid_, const double& mu_, const double& efficiency_, const double& slope_) {
	pid = pid_;
	mu = mu_;
	efficiency = efficiency_;
	slope = slope_;
}

PrimarySource::~PrimarySource() {
}

double PrimarySource::get(const double& E) const {
	double value = 0;
	if (efficiency > 0) {
		double v = beta(E) * c_light;
		double p = beta(E) * pid.get_A() * (E + proton_mass_c2) / c_light;
		value = pid.get_A() * pow2(p) * efficiency * L_SN_surface;
		value /= mu * v * Gamma_Integral(slope) * c_light * pow4(proton_mass_c);
		value *= std::pow(p / proton_mass_c, -slope);
	}
	return value;
}
