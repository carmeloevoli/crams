#include "cgs.h"
#include "primary.h"
#include "utilities.h"

PrimarySource::PrimarySource() {
}

PrimarySource::PrimarySource(const PID& pid, const double& epsilon, const Params& params) {
	A = pid.get_A();
	Z = pid.get_Z();
	slope = (pid == H1) ? params.H_slope : params.nuclei_slope;
	if (epsilon > 0.) {
		double L_SN_surface = cgs::E_SN * cgs::sn_rate / M_PI / pow2(cgs::galaxy_size);
		factor = (double) A * epsilon * L_SN_surface;
		factor /= params.mu * Gamma_Integral(slope) * pow2(cgs::proton_mass_c2);
	}
}

PrimarySource::~PrimarySource() {
	std::cout << "delete primary source for particle " << A << " " << Z << "\n";
}

double PrimarySource::get(const double& T) const {
	double value = 0.;
	if (factor > 0.) {
		double beta = beta_func(T);
		double pc = pc_func(A, T);
		value = factor / beta;
		value *= std::pow(pc / cgs::proton_mass_c2, 2. - slope);
	}
	return value;
}
