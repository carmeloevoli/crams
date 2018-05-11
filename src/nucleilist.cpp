#include "nucleilist.h"

NucleiList::NucleiList() {
}

NucleiList::~NucleiList() {
}

void NucleiList::add_nucleus(const PID& pid, const double& efficiency, const double& gamma) {
	std::pair<double, double> p = std::make_pair(efficiency, gamma);
	list[pid] = p;
}
