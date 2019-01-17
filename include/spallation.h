#ifndef _INCLUDE_SPALLATION_H_
#define _INCLUDE_SPALLATION_H_

#include <map>
#include <string>
#include <vector>
#include "cgs.h"
#include "pid.h"

class SpallationXsecs {
public:
	SpallationXsecs(const PID& fragment);
	virtual ~SpallationXsecs();

	double get(const PID& projectile, const double& T) const;

protected:
	void read_table();

protected:
	PID _fragment;
	std::map<PID, std::vector<double> > _table;
	std::vector<double> _T;
	std::string _table_filename = "xsecs_0.1_100_1.1.txt";
	double _T_min = 0.1 * cgs::GeV;
	size_t _T_size = 100;
	double _T_ratio = 1.1;
};

#endif /* INCLUDE_SPALLATION_H_ */
