#ifndef INCLUDE_PARTICLE_H_
#define INCLUDE_PARTICLE_H_

#include <algorithm>
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

	double get_efficiency() const {
		return _efficiency;
	}

	Grammage* getX() const {
		return _X;
	}

	PID get_pid() const {
		return _pid;
	}

	double get_I(const size_t& i) const {
		return _I_T.at(i);
	}

	inline bool isDone() const {
		return _isDone;
	}

	bool& setDone() {
		return _isDone;
	}

	void build_grammage(const Params& params);
	void build_snr_source(const Params& params);
	void build_inelastic_Xsec(const Params& params);
	void build_losses(const Params& params);
	void build_secondary_source(const std::vector<Particle>& particles, const Params& params);
	void build_tertiary_source(const std::vector<Particle>& particles);
	void build_grammage_at_source(const std::vector<Particle>& particles, const Params& params);
	bool run(const std::vector<double>& T);
	double I_T_interpol(const double& T) const;
	double I_R_LIS(const double& R) const;
	double I_R_TOA(const double& R, const double& modulation_potential) const;
	void dump() const;

public:
	std::string make_filename() const;
	double Q(const double& T);
	double Lambda_1(const double& T);
	double Lambda_2(const double& T);
	double internal_integrand(const double& T_second);
	double ExpIntegral(const double& T, const double& T_prime);
	double external_integrand(const double& T_prime, const double& T);
	double compute_integral(const double& T);

protected:
	bool _isDone = false;
	bool _doGrammageAtSource = false;
	double _efficiency = 0;
	std::vector<double> _T;
	std::vector<double> _I_T;
	PID _pid;
	Grammage* _X = nullptr;
	SnrSource* _Q = nullptr;
	SourceTerm* _Q_sec = nullptr;
	SourceTerm* _Q_ter = nullptr;
	SourceTerm* _Q_Xs = nullptr;
	InelasticXsec* _sigma = nullptr;
	Losses* _dEdx = nullptr;
};

typedef std::vector<Particle> Particles;

struct ptr_Particle {
	bool isPresent;
	std::vector<Particle>::iterator it;
};

#endif /* INCLUDE_PARTICLE_H_ */
