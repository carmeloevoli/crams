#ifndef INCLUDE_PARAMS_H_
#define INCLUDE_PARAMS_H_

#include <string>
#include "mks.h"

static const double E_SN = 1e51 * erg;
static const double sn_rate = 1. / (30. * year);
static const double galaxy_size = 10. * kpc;
static const double L_SN_surface = E_SN * sn_rate / M_PI / pow2(galaxy_size);
static const double He_abundance = 0.15;
static const double mean_ism_mass = proton_mass * (1. + 4 * He_abundance) / (1. + He_abundance);

template<typename T>
class Param {
public:
	Param(const T& v) :
			value(v) {
	}

	virtual ~Param() {
	}

	void set(const T& v) {
		value = v;
	}

	T get() const {
		return value;
	}

protected:
	T value = T();
};

class Params {
public:
	Params() {
	}

	virtual ~Params() {
	}

	void print() {
		//TODO do this
	}

	Param<double> galaxy_H = Param<double>(4. * kpc);
	Param<double> mu = Param<double>(2.4 * mgram / cm2);
	Param<double> D_0 = Param<double>(1e30 * pow2(cm) / sec);
	Param<double> delta_low = Param<double>(0.7);
	Param<double> R_b = Param<double>(500. * GeV);
	Param<double> delta_hi = Param<double>(0.3);
	Param<double> delta_s = Param<double>(0.1);
	Param<double> h_gas = Param<double>(150 * pc);
	Param<double> potential = Param<double>(100 * MeV);
	Param<double> E_min = Param<double>(0.1 * GeV);
	Param<double> E_max = Param<double>(100 * TeV);
	Param<int> E_size = Param<int>(6 * 32);
	double B_0 = muG;
	double ion_number_density = 0.02 / cm3; // TODO make param
    Param<double> v_A = Param<double>(B_0 / std::sqrt(vacuum_permeability * proton_mass * ion_number_density));

	Param<std::string> out_name = Param<std::string>("test");
};

#endif /* INCLUDE_PARAMS_H_ */
