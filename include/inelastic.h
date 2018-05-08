#ifndef INCLUDE_SPALLATION_H_
#define INCLUDE_SPALLATION_H_

#include <cmath>

#include "mks.h"
#include "pid.h"

class InelasticXsec {
public:
	InelasticXsec();
	InelasticXsec(const PID& pid_);
	virtual ~InelasticXsec();
	double get(const double& E) const;

protected:
	PID pid;
	double E_threshold = 0.2797 * GeV;

private:
	double sigma_pp(const double& E) const;
	double sigma_ST(const double& E, const int& A_k) const;
};

#endif /* INCLUDE_SPALLATION_H_ */
