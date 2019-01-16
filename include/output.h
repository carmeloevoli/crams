#ifndef INCLUDE_OUTPUT_H_
#define INCLUDE_OUTPUT_H_

#include <vector>
#include "axis.h"
#include "particle.h"

class OutputManager {
public:
	OutputManager(const LogAxis& T, const std::vector<Particle>& particles);
	virtual ~OutputManager();
	void dump_spectra(double R_min, double R_max, size_t R_size) const;
//	void fill_rigidity(const double& R_min, const double& R_max, const size_t& R_size);
//	void add_nucleus(const PID& pid, const double& efficiency, const double& gamma);
//	void fill_particles(const NucleiList& nucleilist);
//	void dump();
private:
	LogAxis _T;
	std::vector<Particle> _particles;
};

#endif /* INCLUDE_OUTPUT_H_ */
