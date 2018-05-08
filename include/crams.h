#ifndef INCLUDE_CRAMS_H_
#define INCLUDE_CRAMS_H_

#include <vector>

#include "params.h"
#include "particle.h"

class CRAMS {
public:
	CRAMS(const Params& par_);
	virtual ~CRAMS();
	void fill_energy(const double& E_min, const double& E_max, const size_t& E_size);
	void add_nucleus(const PID& pid, const double& efficiency, const double& gamma);
	void fill_particles();
	void test();
	void run();
	void dump();

protected:
	Params par;
	std::vector<double> E;
	std::map<PID, std::pair<double, double> > nucleilist;
	std::vector<Particle> particles; // TODO make it map?
};

#endif /* INCLUDE_CRAMS_H_ */
