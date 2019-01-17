#include "spallation.h"
#include "cgs.h"

SpallationXsecs::SpallationXsecs(const PID& fragment) :
		_fragment(fragment) {
}

SpallationXsecs::~SpallationXsecs() {
}

double SpallationXsecs::get(const PID& projectile, const double& T) const {
	return 1. * cgs::mbarn;
}

void SpallationXsecs::read_table(std::string filename) {
}
