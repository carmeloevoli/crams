#ifndef INCLUDE_LOSSES_H_
#define INCLUDE_LOSSES_H_

#include <cmath>
#include <string>
#include <gsl/gsl_deriv.h>
#include "cgs.h"

#include "params.h"
#include "pid.h"
#include "utilities.h"

class Losses {
public:
	Losses();
	Losses(const PID& pid, const Params& params);
	virtual ~Losses();
	double get(const double& T) const;
	double dE_dx_adiabatic(const double& T) const;
	double dE_dx_ionization(const double& T) const;
	double dE_dt_adiabatic(const double& T) const;
	double dE_dt_ionization(const double& T) const;
	double get_derivative(const double& T);

protected:
	int A = 0;
	int Z = 0;
	double factor_ad = 0;
//	double v_A = 0;
	double mu = 0;
//	double rho_0 = 0;
//	double E_H = 19.e-6 * cgs::MeV;
//	double E_He = 44.e-6 * cgs::MeV;
//	double pi_re2_mc2_c = M_PI * pow2(cgs::electron_radius) * cgs::c_light * cgs::electron_mass_c2;
};

#endif /* INCLUDE_LOSSES_H_ */
