#ifndef INCLUDE_SPALLATION_H_
#define INCLUDE_SPALLATION_H_

#include <cmath>

#include "cgs.h"
#include "params.h"
#include "pid.h"

double sigma_pp(const double& T);

namespace Tripathi99 {
double inelastic_sigma(int A_p, int Z_p, int A_t, int Z_t, double T_n);
} /* namespace Tripathi99 */

class InelasticXsec {
public:
	InelasticXsec();
	InelasticXsec(const PID& pid);
	virtual ~InelasticXsec();
	double get(const double& T) const;

protected:
	int A = 0;
	int Z = 0;
private:
	double sigma_ST(const double& T) const;
};

#endif /* INCLUDE_SPALLATION_H_ */
