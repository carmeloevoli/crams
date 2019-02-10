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
	Chi2(const Particles& particles, const double& phi);
	virtual ~Chi2();
	double compute_chi2(const double& R_min, const double& R_max = 10. * cgs::TeV) const;

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
	void read_datafile(const std::string& filename, const double& units = 1.0);
	virtual double get_model(const double& R, const double& phi) const = 0;

protected:
	double _chi2 = 0;
	double _phi = 0;
	std::vector<data_point> _data;
	Particles _particles;
};

class Chi2_C: public Chi2 {
public:
	Chi2_C(const Particles& particles, const double& phi) :
			Chi2(particles, phi) {
		constexpr double units = 1. / (cgs::GeV * cgs::m2 * cgs::sec * cgs::sr);
		read_datafile("data/C_AMS02_rig.txt", units);
	}
protected:
	double get_model(const double& R, const double& phi) const override;
	ptr_Particle ptr_C12 = find_ptr(C12);
	ptr_Particle ptr_C13 = find_ptr(C13);
	ptr_Particle ptr_C14 = find_ptr(C14);
};

class Chi2_N: public Chi2 {
public:
	Chi2_N(const Particles& particles, const double& phi) :
			Chi2(particles, phi) {
		constexpr double units = 1. / (cgs::GeV * cgs::m2 * cgs::sec * cgs::sr);
		read_datafile("data/N_AMS02_rig.txt", units);
	}
protected:
	double get_model(const double& R, const double& phi) const override;
	ptr_Particle ptr_N14 = find_ptr(N14);
	ptr_Particle ptr_N15 = find_ptr(N15);
};

class Chi2_O: public Chi2 {
public:
	Chi2_O(const Particles& particles, const double& phi) :
			Chi2(particles, phi) {
		constexpr double units = 1. / (cgs::GeV * cgs::m2 * cgs::sec * cgs::sr);
		read_datafile("data/O_AMS02_rig.txt", units);
	}
protected:
	double get_model(const double& R, const double& phi) const override;
	ptr_Particle ptr_O16 = find_ptr(O16);
	ptr_Particle ptr_O17 = find_ptr(O17);
	ptr_Particle ptr_O18 = find_ptr(O18);
};

class Chi2_BC: public Chi2 {
public:
	Chi2_BC(const Particles& particles, const double& phi) :
			Chi2(particles, phi) {
		read_datafile("data/BC_AMS02_rig.txt");
	}
protected:
	double get_model(const double& R, const double& phi) const override;
	ptr_Particle ptr_C12 = find_ptr(C12);
	ptr_Particle ptr_C13 = find_ptr(C13);
	ptr_Particle ptr_C14 = find_ptr(C14);
	ptr_Particle ptr_B10 = find_ptr(B10);
	ptr_Particle ptr_B11 = find_ptr(B11);
};

class Chi2_CO: public Chi2 {
public:
	Chi2_CO(const Particles& particles, const double& phi) :
			Chi2(particles, phi) {
		read_datafile("data/CO_AMS02_rig.txt");
	}
protected:
	double get_model(const double& R, const double& phi) const override;
	ptr_Particle ptr_C12 = find_ptr(C12);
	ptr_Particle ptr_C13 = find_ptr(C13);
	ptr_Particle ptr_C14 = find_ptr(C14);
	ptr_Particle ptr_O16 = find_ptr(O16);
	ptr_Particle ptr_O17 = find_ptr(O17);
	ptr_Particle ptr_O18 = find_ptr(O18);
};

class Chi2_HeO: public Chi2 {
public:
	Chi2_HeO(const Particles& particles, const double& phi) :
			Chi2(particles, phi) {
		read_datafile("data/HeO_AMS02_rig.txt");
	}
protected:
	double get_model(const double& R, const double& phi) const override;
	ptr_Particle ptr_He3 = find_ptr(He3);
	ptr_Particle ptr_He4 = find_ptr(He4);
	ptr_Particle ptr_O16 = find_ptr(O16);
	ptr_Particle ptr_O17 = find_ptr(O17);
	ptr_Particle ptr_O18 = find_ptr(O18);
};

class Chi2_He: public Chi2 {
public:
	Chi2_He(const Particles& particles, const double& phi) :
			Chi2(particles, phi) {
		constexpr double units = 1. / (cgs::GeV * cgs::m2 * cgs::sec * cgs::sr);
		read_datafile("data/He_AMS02_rig.txt", units);
	}
protected:
	double get_model(const double& R, const double& phi) const override;
	ptr_Particle ptr_He3 = find_ptr(He3);
	ptr_Particle ptr_He4 = find_ptr(He4);
};

class Chi2_H: public Chi2 {
public:
	Chi2_H(const Particles& particles, const double& phi) :
			Chi2(particles, phi) {
		constexpr double units = 1. / (cgs::GeV * cgs::m2 * cgs::sec * cgs::sr);
		read_datafile("data/H_AMS02_rig.txt", units);
	}
protected:
	double get_model(const double& R, const double& phi) const override;
	ptr_Particle ptr_H1 = find_ptr(H1);
	ptr_Particle ptr_H2 = find_ptr(H2);
	ptr_Particle ptr_H1_ter = find_ptr(H1_ter);
};

#endif /* INCLUDE_CHI2_H_ */
