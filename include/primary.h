#ifndef INCLUDE_PRIMARY_H_
#define INCLUDE_PRIMARY_H_

#include <cmath>

#include "params.h"
#include "pid.h"
#include "utilities.h"

class PrimarySource {
public:
	PrimarySource();
	PrimarySource(const PID& pid_, const double& mu_, const double& efficiency_, const double& slope_);
	virtual ~PrimarySource();
	double get(const double& E) const;

protected:
	PID pid;
	double mu = 0;
	double efficiency = 0;
	double slope = 0;
};

#endif /* INCLUDE_PRIMARY_H_ */
