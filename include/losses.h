#ifndef INCLUDE_LOSSES_H_
#define INCLUDE_LOSSES_H_

#include <cmath>
#include <string>
#include <gsl/gsl_deriv.h>

#include "mks.h"
#include "params.h"
#include "pid.h"
#include "utilities.h"

class Losses {
public:
	Losses();
	Losses(const PID& pid_, const Params& par);
	virtual ~Losses();
	double adiabatic(const double& E) const;
	double dE_dt_adiabatic(const double& E) const;
	double dE_dt_ionization(const double& E) const;
	double ionization(const double& E) const;
	double dE_dx(const double& E) const;
	double get_derivative(const double& E);

private:

protected:
	PID pid;
	double v_A = 0;
	double mu = 0;
	double rho_0 = 0;
	double E_H = 19.e-6 * MeV;
	double E_He = 44.e-6 * MeV;
	double pi_re2_mc2_c = M_PI * pow2(electron_radius) * c_light * electron_mass_c2;
};

#endif /* INCLUDE_LOSSES_H_ */
