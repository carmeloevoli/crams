#ifndef INCLUDE_SPALLATION_H_
#define INCLUDE_SPALLATION_H_

#include <cmath>

#include "cgs.h"
#include "params.h"
#include "pid.h"

class InelasticXsec {
public:
	InelasticXsec();
	InelasticXsec(const PID& pid);
	virtual ~InelasticXsec();
	double get(const double& T) const;

protected:
	int A = 0;
	int Z = 0;
	double E_threshold = 0.2797 * cgs::GeV;

private:
	double sigma_pp(const double& T) const;
	double sigma_ST(const double& T) const;
};

#endif /* INCLUDE_SPALLATION_H_ */
