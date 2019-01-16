#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "particle.h"
#include "axis.h"
#include "cgs.h"
#include "utilities.h"

#define LIMIT 10000

Particle::Particle() {
}

Particle::Particle(PID pid, double efficiency) :
		_pid(pid), _efficiency(efficiency) {
}

std::string Particle::make_filename() {
	std::string filename = "particle";
	filename += "_" + std::to_string(_pid.get_Z());
	filename += "_" + std::to_string(_pid.get_A()) + ".log";
	return filename;
}

void Particle::dump() {
	auto filename = make_filename();
	std::ofstream outfile(filename);
	outfile << std::scientific;
	for (double T_ = 0.01 * cgs::GeV; T_ < 1e3 * cgs::TeV; T_ *= 1.1) {
		outfile << T_ / cgs::GeV << "\t";
		double E = T_ + cgs::proton_mass_c2;
		double R = pc_func(_pid.get_A(), T_) / (double) _pid.get_Z();
		outfile << R / cgs::GeV << "\t";
		outfile << Q->get(T_) << "\t";
		outfile << X->get(T_) / (cgs::gram / cgs::cm2) << "\t";
		outfile << cgs::mean_ism_mass / sigma->get(T_) / (cgs::gram / cgs::cm2) << "\t";
		outfile << T_ / -dEdx->get(T_) / (cgs::gram / cgs::cm2) << "\t";
		outfile << "\n";
	}
	outfile.close();
}

void Particle::clear() {
	if (X != nullptr)
		delete X;
	if (Q != nullptr)
		delete Q;
	if (Q_sec != nullptr)
		delete Q_sec;
	if (sigma != nullptr)
		delete sigma;
	if (dEdx != nullptr)
		delete dEdx;
}

Particle::~Particle() {
}

double Particle::get_I_T_interpol(const double& T) const {
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	size_t size = _T.size();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, size);
	double logT[size], logI_T[size];
	for (size_t i = 0; i < size; ++i) { // TODO make more efficient
		logT[i] = std::log(_T[i]);
		logI_T[i] = (_I_T[i] > 0.) ? std::log(_I_T[i]) : -30;
	}
	gsl_spline_init(spline, logT, logI_T, size);
	double value = gsl_spline_eval(spline, std::log(T), acc);
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	return std::exp(value);
}

double Particle::get_I_R_LIS(const double& R) const {
	double value = 0;
	constexpr double mp_2 = pow2(cgs::proton_mass_c2);
	double E_2 = pow2(R * _pid.get_Z_over_A()) + mp_2;
	double T_ = std::sqrt(E_2) - cgs::proton_mass_c2;
	if (T_ > _T.front() && T_ < _T.back()) {
		double Z_A_squared = pow2(_pid.get_Z_over_A());
		double dTdR = R * Z_A_squared;
		dTdR /= std::sqrt(Z_A_squared * pow2(R) + mp_2);
		value = get_I_T_interpol(T_) * dTdR;
	}
	return value;
}

double Particle::get_I_R_TOA(const double& R, const double& modulation_potential) const {
	double value = 0;
	constexpr double mp_2 = pow2(cgs::proton_mass_c2);
	double E_2 = pow2(R * _pid.get_Z_over_A()) + pow2(cgs::proton_mass_c2);
	double T_now = std::sqrt(E_2) - cgs::proton_mass_c2;
	if (T_now > _T.front() && T_now < _T.back()) {
		double Phi = _pid.get_Z_over_A() * modulation_potential;
		double T_ISM = std::min(T_now + Phi, _T.back());
		double Z_A_squared = pow2(_pid.get_Z_over_A());
		double dTdR = R * Z_A_squared;
		dTdR /= std::sqrt(Z_A_squared * pow2(R) + mp_2);
		double factor = (pow2(T_now) - mp_2) / (pow2(T_ISM) - mp_2);
		value = factor * get_I_T_interpol(T_ISM) * dTdR;
	}
	return value;
}

int Particle::run_spectrum(const LogAxis& T) {
	_T = T.get_axis();
	for (auto& T_now : _T) {
		_I_T.push_back(compute_integral(T_now));
	}
	assert(_I_T.size() == _T.size());
	return 0;
}

double Particle::Lambda_1(const double& T) {
	double value = 1. / X->get(T) + sigma->get(T) / cgs::mean_ism_mass + dEdx->get_derivative(T);
	return value;
}

double Particle::Lambda_2(const double& T) {
	double value = dEdx->get(T);
	return std::fabs(value);
}

double compute_integral_qags(gsl_integration_workspace * w, gsl_function * F, double x_lo,
		double x_hi) {
	double result, error;
	gsl_integration_qags(F, x_lo, x_hi, 0, 1e-5, LIMIT, w, &result, &error);
	return result;
}

double Particle::internal_integrand(const double& T_second) {
	return Lambda_1(T_second) / Lambda_2(T_second);
}

double gslParticleClassExpWrapper(double x, void * pp) {
	double E_second = std::exp(x);
	Particle * particle = (Particle *) pp;
	return E_second * particle->internal_integrand(E_second);
}

double Particle::ExpIntegral(const double& T, const double& T_prime) {
	double result = 0;
	gsl_function F;
	F.params = this;
	F.function = &gslParticleClassExpWrapper;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(LIMIT);
	result = compute_integral_qags(w, &F, std::log(T), std::log(T_prime));
	gsl_integration_workspace_free(w);
	return result;
}

double Particle::external_integrand(const double& T_prime, const double& T) {
	return Q->get(T_prime) / Lambda_2(T_prime) * std::exp(-ExpIntegral(T, T_prime));
}

struct gsl_f_pars {
	double T;
	Particle * pt_Particle;
};

double gslParticleClassWrapper(double x, void * pp) {
	gsl_f_pars *p = (gsl_f_pars *) pp;
	double T = p->T;
	double T_prime = std::exp(x);
	return T_prime * p->pt_Particle->external_integrand(T_prime, T);
}

double Particle::compute_integral(const double& T) {
	double result = 0;
	gsl_f_pars pars = { T, this };
	gsl_function F;
	F.params = &pars;
	F.function = &gslParticleClassWrapper;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(LIMIT);
	result = compute_integral_qags(w, &F, std::log(T), std::log(1e3 * T));
	gsl_integration_workspace_free(w);
	return result;
}

//void Particle::dump_timescales() const {
//	std::cout << std::scientific;
//	for (double E = 0.1 * GeV; E < 10 * PeV; E *= 1.1) {
//		std::cout << E / GeV << "\t";
//		std::cout << X.get(E) / (gram / cm2) << "\t";
//		std::cout << X.diffusion_escape_time(E) / year << "\t";
//		std::cout << X.advection_escape_time() / year << "\t";
//		std::cout << std::pow(E, 2.2) * Q.get(E) << "\t";
//		std::cout << -(E / b.dE_dt_adiabatic(E)) / year << "\t";
//		std::cout << -(E / b.dE_dt_ionization(E)) / year << "\t";
//		std::cout << sigma_in.get(E) / mbarn << "\t";
//		std::cout << "\n";
//	}
//}

#undef LIMIT
