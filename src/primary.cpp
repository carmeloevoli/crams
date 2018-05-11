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

double PrimarySource::get(const double& T) const {
	double value = 0;
	if (efficiency > 0) {
		double p_c = beta(T) * pid.get_A() * (T + proton_mass_c2);
		value = pid.get_A() * pow2(p_c) * efficiency * L_SN_surface;
		value /= mu * beta(T) * Gamma_Integral(slope) * pow4(proton_mass_c2);
		value *= std::pow(p_c / proton_mass_c2, -slope);
	}
	return value;
}
