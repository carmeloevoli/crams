#ifndef INCLUDE_GRAMMAGE_H_
#define INCLUDE_GRAMMAGE_H_

#include <cmath>
#include "pid.h"
#include "params.h"
#include "utilities.h"

class Grammage {
public:
	Grammage();
	Grammage(const PID& pid_, const Params& par);
	virtual ~Grammage();
	double D(const double& beta_, const double& E) const;
	double diffusion_escape_time(const double& E) const;
	double advection_escape_time() const;
	double get(const double& E) const;

protected:
	PID pid;
	double mu = 0;
	double v_A = 0;
	double H = 0;
	double R_b = 0;
	double delta_hi = 0;
	double delta_low = 0;
	double delta_s = 0;
	double D_0 = 0;
	double rho_0 = 0;
};

#endif /* INCLUDE_GRAMMAGE_H_ */
