#include "secondary.h"

SourceTerm::SourceTerm() {
}

SourceTerm::SourceTerm(const std::vector<double>& T, const std::vector<double>& Q) :
		_T(T), _Q(Q) {
}

SourceTerm::~SourceTerm() {
}

double SourceTerm::get(const double& T) const {
	return LinearInterpolatorLog(_T, _Q, T);
}
