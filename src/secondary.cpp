#include "secondary.h"

SecondarySource::SecondarySource() {
}

SecondarySource::SecondarySource(const std::vector<double>& T, const std::vector<double>& Q) :
		_T(T), _Q(Q) {
}

SecondarySource::~SecondarySource() {
}

double SecondarySource::get(const double& T) const {
	return LinearInterpolatorLog(_T, _Q, T);
}
