#ifndef INCLUDE_PRIMARY_H_
#define INCLUDE_PRIMARY_H_

#include "params.h"
#include "pid.h"

class SnrSource {
public:
	SnrSource();
	SnrSource(const PID& pid, const double& epsilon, const Params& params);
	virtual ~SnrSource();
	double get(const double& T) const;

protected:
	int A = 0;
	int Z = 0;
	double slope = 0.;
	double factor = 0.;
};

#endif /* INCLUDE_PRIMARY_H_ */
