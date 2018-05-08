#ifndef INCLUDE_PARTICLE_H_
#define INCLUDE_PARTICLE_H_

#include <cmath>
#include <iomanip>
#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include "grammage.h"
#include "inelastic.h"
#include "losses.h"
#include "mks.h"
#include "PID.h"
#include "primary.h"
#include "utilities.h"

class Particle {
public:
	Particle();
	Particle(PID pid_, double efficiency_, double slope_, Params par);
	virtual ~Particle();
	void dump_inputs() const;
	double get_I(const size_t& i) const;
	double get_modulated(const size_t& i) const;
	PID get_pid() const;
	void compute_spectrum(const std::vector<double>& Ek);
	void modulate(const std::vector<double>& Ek);
	double external_integrand(const double& E_prime, const double& E);
	double internal_integrand(const double& E_second);
	double ExpIntegral(const double& E, const double& E_prime);

private:
	double compute_integral(const double& E);
	double Lambda_1(const double& E);
	double Lambda_2(const double& E);

protected:
	std::vector<double> I_Ek;
	std::vector<double> I_Ek_mod;
	PID pid;
	Grammage X;
	InelasticXsec sigma_in;
	PrimarySource Q;
	Losses b;
	double potential = 0;
};

#endif /* INCLUDE_PARTICLE_H_ */
