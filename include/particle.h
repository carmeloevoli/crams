#ifndef INCLUDE_PARTICLE_H_
#define INCLUDE_PARTICLE_H_

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "axis.h"
#include "cgs.h"
#include "grammage.h"
#include "inelastic.h"
#include "losses.h"
#include "pid.h"
#include "primary.h"

class Particle {
public:
	Particle();
	Particle(PID pid_, double efficiency_);
	virtual ~Particle();
	void clear();

	bool operator==(const Particle &other) const {
		return pid == other.pid;
	}

	//void dump_timescales() const;
	//void compute_spectrum(const std::vector<double>& T);
	//void modulate(const std::vector<double>& T, const std::vector<double>& R);
	//double external_integrand(const double& E_prime, const double& E);
	//double internal_integrand(const double& E_second);
	//double ExpIntegral(const double& E, const double& E_prime);

	const std::vector<double>& get_I() const {
		return I_T;
	}

	PID get_pid() const {
		return pid;
	}

	double get_I(const size_t& i) const {
		return I_T.at(i);
	}

	void build_grammage(const Params& params) {
		X = new Grammage(pid, params);
	}

	void build_primary_source(const Params& params) {
		Q = new PrimarySource(pid, efficiency, params);
	}

	void build_inelastic_Xsec(const Params& params) {
		sigma = new InelasticXsec(pid);
	}

	void build_losses(const Params& params) {
		dEdx = new Losses(pid, params);
	}

	double get_X(const double& T_) const {
		return X->get(T_);
	}

	double get_Q(const double& T_) const {
		return Q->get(T_);
	}

	double get_sigma_m(const double& T_) const {
		return sigma->get(T_);
	}

	double get_dEdx(const double& T_) const {
		return dEdx->get(T_);
	}

	double get_I_T_interpol(const double& T_) const;
	double get_I_R_LIS(const double& R) const;
	double get_I_R_TOA(const double& R, const double& modulation_potential) const;
	double f(double T_, double Y);
	int run(const LogAxis& T_);
	int run_gsl(const LogAxis& T_);
	int run_spectrum(const LogAxis& T_);
	void dump();

public:
	std::string make_filename();
	double Lambda_1(const double& T);
	double Lambda_2(const double& T);
	double internal_integrand(const double& T_second);
	double ExpIntegral(const double& T, const double& T_prime);
	double external_integrand(const double& T_prime, const double& T);
	double compute_integral(const double& T);

protected:
	std::vector<double> T;
	std::vector<double> I_T;
	PID pid;
	double efficiency = 0;
	Grammage* X = nullptr;
	PrimarySource* Q = nullptr;
	InelasticXsec* sigma = nullptr;
	Losses* dEdx = nullptr;
};

#endif /* INCLUDE_PARTICLE_H_ */
