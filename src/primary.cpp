#include "cgs.h"
#include "primary.h"
#include "utilities.h"

SnrSource::SnrSource() {
}

SnrSource::~SnrSource() {
#ifdef DEBUG
	std::cout << "delete primary source for particle " << A << " " << Z << "\n";
#endif
}

SnrSource::SnrSource(const PID& pid, const double& epsilon, const Params& params) {
	A = pid.get_A();
	Z = pid.get_Z();
	if (pid.get_Z() == 1)
		slope = params.H_slope;
	else if (pid.get_Z() == 2)
		slope = params.He_slope;
	else
		slope = params.nuclei_slope;
	if (epsilon > 0.) {
		double L_SN_surface = cgs::E_SN * cgs::sn_rate / M_PI / pow2(cgs::galaxy_size);
		factor = (double) A * epsilon * L_SN_surface;
		factor /= params.mu * Gamma_Integral(slope) * pow2(cgs::proton_mass_c2);
	}
}

double SnrSource::get(const double& T) const {
	double value = 0.;
	if (factor > 0.) {
		double beta = beta_func(T);
		double pc = pc_func(A, T);
		value = factor / beta;
		value *= std::pow(pc / cgs::proton_mass_c2, 2. - slope);
	}
	return value;
}
