#ifndef INCLUDE_PRIMARY_H_
#define INCLUDE_PRIMARY_H_

#include "params.h"
#include "pid.h"

class PrimarySource {
public:
	PrimarySource();
	PrimarySource(const PID& pid, const double& epsilon, const Params& params);
	virtual ~PrimarySource();
	double get(const double& T) const;

protected:
	int A = 0;
	int Z = 0;
	double slope = 0.;
	double factor = 0.;
};

#endif /* INCLUDE_PRIMARY_H_ */
