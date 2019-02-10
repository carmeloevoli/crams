#ifndef INCLUDE_SECONDARY_H_
#define INCLUDE_SECONDARY_H_

#include <vector>
#include "utilities.h"

class SourceTerm {
public:
	SourceTerm();
	SourceTerm(const std::vector<double>& T, const std::vector<double>& Q);
	virtual ~SourceTerm();
	double get(const double& T) const;

private:
	std::vector<double> _T;
	std::vector<double> _Q;
};

#endif /* INCLUDE_SECONDARY_H_ */
