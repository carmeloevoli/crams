#ifndef INCLUDE_OUTPUT_H_
#define INCLUDE_OUTPUT_H_

#include <vector>
#include "particle.h"

class OutputManager {
public:
	OutputManager(const Particles& particles, const double phi);
	virtual ~OutputManager();
	void dump_spectra(double R_min, double R_max, size_t R_size) const;
	void dump_ratio(double R_min, double R_max, size_t R_size) const;
	void dump_heavy_spectra(double R_min, double R_max, size_t R_size) const;
private:
	Particles _particles;
	double _phi = 0;

	ptr_Particle ptr_H1 = find_ptr(H1);
	ptr_Particle ptr_B10 = find_ptr(B10);
	ptr_Particle ptr_B11 = find_ptr(B11);
	ptr_Particle ptr_C12 = find_ptr(C12);
	ptr_Particle ptr_C13 = find_ptr(C13);
	ptr_Particle ptr_C14 = find_ptr(C14);
	ptr_Particle ptr_N14 = find_ptr(N14);
	ptr_Particle ptr_N15 = find_ptr(N15);
	ptr_Particle ptr_O16 = find_ptr(O16);
	ptr_Particle ptr_O17 = find_ptr(O17);
	ptr_Particle ptr_O18 = find_ptr(O18);
	ptr_Particle find_ptr(const PID& pid);

	double H(const double& R) const;
	double B(const double& R) const;
	double C(const double& R) const;
	double N(const double& R) const;
	double O(const double& R) const;
	double C12_C13(const double& R) const;
};

#endif /* INCLUDE_OUTPUT_H_ */
