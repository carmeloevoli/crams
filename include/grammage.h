#ifndef INCLUDE_GRAMMAGE_H_
#define INCLUDE_GRAMMAGE_H_

#include "pid.h"
#include "params.h"

class Grammage {
public:
	Grammage();
	Grammage(const PID& pid, const Params& params);
	virtual ~Grammage();
	double D(const double& T) const;
	double get(const double& T) const;

protected:
	int A = 0;
	int Z = 0;
	double factor = 0;
	double v_A = 0;
	double H = 0;
	double D_0 = 0;
	double R_b = 0;
	double delta = 0;
	double ddelta = 0;
	double s = 0;
};

#endif /* INCLUDE_GRAMMAGE_H_ */
