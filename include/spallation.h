#ifndef _INCLUDE_SPALLATION_H_
#define _INCLUDE_SPALLATION_H_

#include <map>
#include <string>
#include <vector>
#include "pid.h"

class SpallationXsecs {
public:
	SpallationXsecs(const PID& fragment);
	virtual ~SpallationXsecs();

	double get(const PID& projectile, const double& T) const;

protected:
	void read_table(std::string filename);

protected:
	PID _fragment;
	std::map<PID, std::vector<double> > table;
};

#endif /* INCLUDE_SPALLATION_H_ */
