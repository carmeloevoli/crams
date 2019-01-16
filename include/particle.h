#ifndef INCLUDE_PARTICLE_H_
#define INCLUDE_PARTICLE_H_

#include <vector>

#include "axis.h"
#include "grammage.h"
#include "inelastic.h"
#include "losses.h"
#include "pid.h"
#include "primary.h"
#include "secondary.h"

class Particle {
public:
	Particle();
	Particle(PID pid_, double efficiency_);
	virtual ~Particle();
	void clear();

	bool operator==(const Particle &other) const {
		return _pid == other._pid;
	}

	const std::vector<double>& get_I() const {
		return _I_T;
	}

	PID get_pid() const {
		return _pid;
	}

	double get_I(const size_t& i) const {
		return _I_T.at(i);
	}

	void build_grammage(const Params& params) {
		X = new Grammage(_pid, params);
	}

	void build_primary_source(const Params& params) {
		Q = new PrimarySource(_pid, _efficiency, params);
	}

	void build_secondary_source(const Params& params) {
		Q_sec = new SecondarySource();
	}

	void build_inelastic_Xsec(const Params& params) {
		sigma = new InelasticXsec(_pid);
	}

	void build_losses(const Params& params) {
		dEdx = new Losses(_pid, params);
	}

	double get_X(const double& T) const {
		return X->get(T);
	}

	double get_Q(const double& T) const {
		return Q->get(T) + Q_sec->get(T);
	}

	double get_sigma_m(const double& T) const {
		return sigma->get(T);
	}

	double get_dEdx(const double& T) const {
		return dEdx->get(T);
	}

	double get_I_T_interpol(const double& T) const;
	double get_I_R_LIS(const double& R) const;
	double get_I_R_TOA(const double& R, const double& modulation_potential) const;
	void dump();

	//odeint.cpp
	double f(double T, double Y);
	int run(const LogAxis& T);
	int run_gsl(const LogAxis& T);
	int run_spectrum(const LogAxis& T);

public:
	std::string make_filename();
	double Lambda_1(const double& T);
	double Lambda_2(const double& T);
	double internal_integrand(const double& T_second);
	double ExpIntegral(const double& T, const double& T_prime);
	double external_integrand(const double& T_prime, const double& T);
	double compute_integral(const double& T);

protected:
	std::vector<double> _T;
	std::vector<double> _I_T;
	PID _pid;
	double _efficiency = 0;
	Grammage* X = nullptr;
	PrimarySource* Q = nullptr;
	SecondarySource* Q_sec = nullptr;
	InelasticXsec* sigma = nullptr;
	Losses* dEdx = nullptr;
};

#endif /* INCLUDE_PARTICLE_H_ */
