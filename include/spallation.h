#ifndef INCLUDE_SPALLATION_H_
#define INCLUDE_SPALLATION_H_

#include "pid.h"

class SpallationXsecs {
public:
	SpallationXsecs() {
	}

	virtual ~SpallationXsecs() {
	}

	double get(const Channel& channel, const double& T_) {
		return 0;
	}

};

#endif /* INCLUDE_SPALLATION_H_ */
