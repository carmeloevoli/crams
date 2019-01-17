#ifndef INCLUDE_PARTICLE_H_
#define INCLUDE_PARTICLE_H_

#include <vector>

#include "grammage.h"
#include "inelastic.h"
#include "losses.h"
#include "pid.h"
#include "primary.h"
#include "spallation.h"
#include "secondary.h"
#include "utilities.h"

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
		_Q = new PrimarySource(_pid, _efficiency, params);
	}

	inline bool isDone() const {
		return _isDone;
	}

	bool& setDone() {
		return _isDone;
	}

	void build_secondary_source(const std::vector<Particle>& particles) {
		auto xsecs = SpallationXsecs(_pid);
		const size_t size = 100;
		auto T_s = LogAxis(0.1 * cgs::GeV, 10. * cgs::TeV, size);
		std::vector<double> Q_s;
		for (auto& T : T_s) {
			double value = 0;
			for (auto& particle : particles) {
				if (particle.get_pid().get_A() > _pid.get_A() && particle.isDone()) {
					value += xsecs.get(particle.get_pid(), T) * particle.get_I_T_interpol(T);
				}
			}
			value /= cgs::mean_ism_mass;
			Q_s.push_back(value);
		}
		_Q_sec = new SecondarySource(T_s, Q_s);
	}

	void build_inelastic_Xsec() {
		_sigma = new InelasticXsec(_pid);
	}

	void build_losses(const Params& params) {
		_dEdx = new Losses(_pid, params);
	}

	bool run(const std::vector<double>& T);
	double get_I_T_interpol(const double& T) const;
	double get_I_R_LIS(const double& R) const;
	double get_I_R_TOA(const double& R, const double& modulation_potential) const;
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
	bool _isDone = false;
	double _efficiency = 0;
	std::vector<double> _T;
	std::vector<double> _I_T;
	PID _pid;
	Grammage* X = nullptr;
	PrimarySource* _Q = nullptr;
	SecondarySource* _Q_sec = nullptr;
	InelasticXsec* _sigma = nullptr;
	Losses* _dEdx = nullptr;
};

#endif /* INCLUDE_PARTICLE_H_ */
