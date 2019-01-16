#ifndef INCLUDE_PRIMARY_H_
#define INCLUDE_PRIMARY_H_

#include <cmath>

#include "params.h"
#include "pid.h"
#include "utilities.h"

class PrimarySource {
public:
	PrimarySource();
	PrimarySource(const PID& pid, const double& epsilon, const Params& params);
	virtual ~PrimarySource();
	double get(const double& T) const;

protected:
	const double L_SN_surface = cgs::E_SN * cgs::sn_rate / M_PI / pow2(cgs::galaxy_size);
	int A = 0;
	int Z = 0;
	double slope = 0.;
	double factor = 0.;
};

#endif /* INCLUDE_PRIMARY_H_ */
