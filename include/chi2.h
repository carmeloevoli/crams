#ifndef INCLUDE_CHI2_H_
#define INCLUDE_CHI2_H_

#include <string>
#include <vector>
#include "particle.h"
#include "pid.h"

struct data_point {
	double R;
	double F;
	double F_err_low;
	double F_err_high;
};

class Chi2 {
public:
	Chi2();
	Chi2(const Particles& particles, const double& phi, const std::string& filename);
	virtual ~Chi2();
	double compute_chi2(const double& R_min);

	ptr_Particle find_ptr(const PID& pid) {
		auto it = find(_particles.begin(), _particles.end(), Particle(pid, 0));
		bool isPresent = !(it == _particles.end());
		return {isPresent, it};
	}

	void set_phi(const double& modulation_potential) {
		_phi = modulation_potential;
	}

	inline double getPhi() const {
		return _phi;
	}

	inline double getChi2() const {
		return _chi2;
	}

protected:
	void read_datafile(const std::string& filename);
	virtual double get_I_R_TOA(const double& R, const double& phi) = 0;

protected:
	double _chi2 = 0;
	double _phi = 0;
	std::vector<data_point> _data;
	Particles _particles;
};

class Chi2_C: public Chi2 {
public:
	Chi2_C(const Particles& particles, const double& phi, const std::string& filename) :
			Chi2(particles, phi, filename) {
	}
protected:
	double get_I_R_TOA(const double& R, const double& phi) override;
	ptr_Particle ptr_C12 = find_ptr(C12);
	ptr_Particle ptr_C13 = find_ptr(C13);
	ptr_Particle ptr_C14 = find_ptr(C14);
};

class Chi2_O: public Chi2 {
public:
	Chi2_O(const Particles& particles, const double& phi, const std::string& filename) :
			Chi2(particles, phi, filename) {
	}
protected:
	double get_I_R_TOA(const double& R, const double& phi) override;
	ptr_Particle ptr_O16 = find_ptr(O16);
	ptr_Particle ptr_O17 = find_ptr(O17);
	ptr_Particle ptr_O18 = find_ptr(O18);
};

#endif /* INCLUDE_CHI2_H_ */
