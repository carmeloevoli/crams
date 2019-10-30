#ifndef INCLUDE_PARAMS_H_
#define INCLUDE_PARAMS_H_

#include <map>

#include "cgs.h"
#include "pid.h"

typedef std::map<PID, double> List;

class ParticleList {
private:
	List _list;

public:
	ParticleList() {
		//insert(H1_ter, 0.);
		//insert(H1, 7e-2);
		//insert(H2, 0.);
		//insert(He3, 0.);
		//insert(He4, 2.5e-2);
		//insert(Li6, 0.);
		//insert(Li7, 0.);
		insert(Be7, 0);
		insert(Be9, 0);
		insert(B10, 0.);
		insert(Be10, 0);
		insert(B11, 0.);
		insert(C12, 5e-3);
		insert(C13, 0.);
		insert(C14, 0.);
		insert(N14, 1e-4);
		insert(N15, 0.);
		insert(O16, 8e-3);
		insert(O17, 0.);
		insert(O18, 0.);
		insert(Ne20, 1e-3);
		insert(Mg24, 2e-3);
		insert(Si28, 2e-3);
		insert(Fe56, 8e-3);
	}

	const List& get_list() {
		const List& l = _list;
		return l;
	}

	virtual ~ParticleList();
	bool insert(const PID& key, const double& value);
	void set_abundance(const PID& key, const double& value);
	void print();
	void set_from_file(const std::string& filename);
};

class Params {
private:
	double _H = 4. * cgs::kpc;
	double _T_min = 0.1 * cgs::GeV;
	double _T_max = 10. * cgs::TeV;
	double _mu = 2.7 * cgs::mgram / cgs::cm2;
	double _v_A = 10. * cgs::km / cgs::sec;
	double _R_b = 312. * cgs::GeV;
	double _delta = 0.62;
	double _ddelta = 0.15;
	double _smoothness = 0.1;
	double _D_0 = 1.8e28 * cgs::cm2 / cgs::sec;
	double _X_s = 0. * cgs::gram / cgs::cm2;
	double _H_slope = 4.25;
	double _He_slope = 4.25;
	double _nuclei_slope = 4.25;
	double _modulation_potential = 0.7 * cgs::GeV;
	double _xsecs_norm = 1;
	size_t _T_size = 100;
	size_t _id = 0;

public:
	Params() {
	}

	void set_H(const double& _H) {
		this->_H = _H;
	}

	void set_D0(const double& _D_0) {
		this->_D_0 = _D_0;
	}

	void set_vA(const double& _v_A) {
		this->_v_A = _v_A;
	}

	void set_delta(const double& _delta) {
		this->_delta = _delta;
	}

	virtual ~Params();
	void print();
	void set_from_file(const std::string& filename);
	void set_params(const std::string& key, const double& value);

	const double& T_min = _T_min;
	const double& T_max = _T_max;
	const double& H = _H;
	const double& mu = _mu;
	const double& v_A = _v_A;
	const double& R_b = _R_b;
	const double& delta = _delta;
	const double& ddelta = _ddelta;
	const double& smoothness = _smoothness;
	const double& D_0 = _D_0;
	const double& X_s = _X_s;
	const double& H_slope = _H_slope;
	const double& He_slope = _He_slope;
	const double& nuclei_slope = _nuclei_slope;
	const double& modulation_potential = _modulation_potential;
	const double& xsecs_norm = _xsecs_norm;
	const size_t& T_size = _T_size;
	const size_t& id = _id;
};

#endif /* INCLUDE_PARAMS_H_ */
