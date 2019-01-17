#ifndef INCLUDE_OUTPUT_H_
#define INCLUDE_OUTPUT_H_

#include <vector>
#include "particle.h"

class OutputManager {
public:
	OutputManager(const std::vector<Particle>& particles, const double phi);
	virtual ~OutputManager();
	void dump_spectra(double R_min, double R_max, size_t R_size) const;
	void dump_heavy_spectra(double R_min, double R_max, size_t R_size) const;
	void dump_ratio(double R_min, double R_max, size_t R_size) const;
private:
	std::vector<Particle> _particles;
	double _phi = 0;
};

#endif /* INCLUDE_OUTPUT_H_ */
