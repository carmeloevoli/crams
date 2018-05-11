#ifndef INCLUDE_CRAMS_H_
#define INCLUDE_CRAMS_H_

#include <vector>

#include "nucleilist.h"
#include "params.h"
#include "particle.h"

class CRAMS {
public:
	CRAMS(const Params& par_);
	virtual ~CRAMS();
	void fill_energy(const double& E_min, const double& E_max, const size_t& E_size);
	void fill_rigidity(const double& R_min, const double& R_max, const size_t& R_size);
	void add_nucleus(const PID& pid, const double& efficiency, const double& gamma);
	void fill_particles(const NucleiList& nucleilist);
	void run();
	void dump();

	const std::vector<double>& get_T() const {
		return T;
	}

	const std::vector<double>& get_spectrum(const PID& pid) const {
		return particles.at(pid).get_I();
//		if (particles.find(pid) == particles.end()) {
//			std::cout << "particle not found in chi2!\n";
	}

protected:
	Params par;
	std::vector<double> T;
	std::map<PID, Particle> particles; // TODO make it map?
};

#endif /* INCLUDE_CRAMS_H_ */
