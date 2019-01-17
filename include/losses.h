#ifndef INCLUDE_LOSSES_H_
#define INCLUDE_LOSSES_H_

#include "params.h"
#include "pid.h"

class Losses {
public:
	Losses();
	Losses(const PID& pid, const Params& params);
	virtual ~Losses();
	double get(const double& T) const;
	double dE_dx_adiabatic(const double& T) const;
	double dE_dx_ionization(const double& T) const;
	//double dE_dt_adiabatic(const double& T) const;
	//double dE_dt_ionization(const double& T) const;
	double get_derivative(const double& T);

protected:
	int A = 0;
	int Z = 0;
	double factor_ad = 0;
	double mu = 0;
};

#endif /* INCLUDE_LOSSES_H_ */
