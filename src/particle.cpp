#include "particle.h"
#include "axis.h"
#include "utilities.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#define LIMIT 10000

Particle::Particle() {
}

Particle::Particle(PID pid_, double efficiency_) :
		pid(pid_), efficiency(efficiency_) {
}

std::string Particle::make_filename() {
	std::string filename = "particle";
	filename += "_" + std::to_string(pid.get_Z());
	filename += "_" + std::to_string(pid.get_A()) + ".log";
	return filename;
}

void Particle::dump() {
	auto filename = make_filename();
	std::ofstream outfile(filename);
	outfile << std::scientific;
	for (double T_ = 0.01 * cgs::GeV; T_ < 1e3 * cgs::TeV; T_ *= 1.1) {
		outfile << T_ / cgs::GeV << "\t";
		double E = T_ + cgs::proton_mass_c2;
		double R = pc_func(pid.get_A(), T_) / (double) pid.get_Z();
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
	if (sigma != nullptr)
		delete sigma;
	if (dEdx != nullptr)
		delete dEdx;
}

Particle::~Particle() {
}

double Particle::get_I_T_interpol(const double& T_) const {
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	size_t size = T.size();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, size);
	double logT[size], logI_T[size];
	for (size_t i = 0; i < size; ++i) { // TODO make more efficient
		logT[i] = std::log(T[i]);
		logI_T[i] = (I_T[i] > 0.) ? std::log(I_T[i]) : -30;
	}
	gsl_spline_init(spline, logT, logI_T, size);
	double value = gsl_spline_eval(spline, std::log(T_), acc);
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
	return std::exp(value);
}

double Particle::get_I_R_LIS(const double& R) const {
	double value = 0;
	constexpr double mp_2 = pow2(cgs::proton_mass_c2);
	double E_2 = pow2(R * pid.get_Z_over_A()) + mp_2;
	double T_ = std::sqrt(E_2) - cgs::proton_mass_c2;
	if (T_ > T.front() && T_ < T.back()) {
		double Z_A_squared = pow2(pid.get_Z_over_A());
		double dTdR = R * Z_A_squared;
		dTdR /= std::sqrt(Z_A_squared * pow2(R) + mp_2);
		value = get_I_T_interpol(T_) * dTdR;
	}
	return value;
}

double Particle::get_I_R_TOA(const double& R, const double& modulation_potential) const {
	double value = 0;
	constexpr double mp_2 = pow2(cgs::proton_mass_c2);
	double E_2 = pow2(R * pid.get_Z_over_A()) + pow2(cgs::proton_mass_c2);
	double T_ = std::sqrt(E_2) - cgs::proton_mass_c2;
	if (T_ > T.front() && T_ < T.back()) {
		double Phi = pid.get_Z_over_A() * modulation_potential;
		double T_ISM = std::min(T_ + Phi, T.back());
		double Z_A_squared = pow2(pid.get_Z_over_A());
		double dTdR = R * Z_A_squared;
		dTdR /= std::sqrt(Z_A_squared * pow2(R) + mp_2);
		double factor = (pow2(T_) - mp_2) / (pow2(T_ISM) - mp_2);
		value = factor * get_I_T_interpol(T_ISM) * dTdR;
	}
	return value;
}

double Particle::f(double T_, double Y) {
	return Q->get(T_) - Y / -dEdx->get(T_) * (1. / X->get(T_) + sigma->get(T_));
}

int Particle::run(const LogAxis& T_) {
	T = T_.get_axis();
	I_T.resize(T.size());
	double Y_n = 0;
	for (size_t i = T.size() - 1; i > 1; --i) {
		for (double T_n = T.at(i); T_n > T.at(i - 1); T_n /= 1.01) {
			double h = T_n * (1. - 1. / 1.01);
			double k_1 = h * f(T_n, Y_n);
			double k_2 = h * f(T_n - .5 * h, Y_n + .5 * k_1);
			double k_3 = h * f(T_n - .5 * h, Y_n + .5 * k_2);
			double k_4 = h * f(T_n - h, Y_n + k_3);
			Y_n = Y_n + 1. / 6. * (k_1 + 2. * k_2 + 2. * k_3 + k_4);
			std::cout << std::scientific << T_n / cgs::GeV << " " << h / cgs::GeV << " "
					<< Y_n / -dEdx->get(T_n) << "\n";
		}
	}
	return 0;
}

struct gsl_f_params {
	Particle * pt_Particle;
};

int func(double T, const double Y[], double f[], void *params) {
	auto p = *(gsl_f_params *) params;
	Particle * particle = p.pt_Particle;
	T = std::fabs(T);
	double Q = particle->get_Q(T);
	double dEdx = std::fabs(particle->get_dEdx(T));
	double X = particle->get_X(T);
	double sigma_m = particle->get_sigma_m(T);
	f[0] = Q - Y[0] / dEdx * (1. / X + sigma_m);
	return GSL_SUCCESS;
}

int jac(double T, const double Y[], double *dfdY, double dfdT[], void *params) {
	auto p = *(gsl_f_params *) params;
	Particle * particle = p.pt_Particle;
	T = std::fabs(T);
	double dEdx = std::fabs(particle->get_dEdx(T));
	double X = particle->get_X(T);
	double sigma_m = particle->get_sigma_m(T);
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdY, 1, 1);
	gsl_matrix * m = &dfdy_mat.matrix;
	gsl_matrix_set(m, 0, 0, -1. / dEdx * (1. / X + sigma_m));
	dfdT[0] = 0.0;
	return GSL_SUCCESS;
}

int Particle::run_gsl(const LogAxis& T_) {
	T = T_.get_axis();
	gsl_f_params params = { this };
	gsl_odeiv2_system sys = { func, jac, 1, &params };
	double Y[1] = { 0.0 };
	I_T.reserve(T.size());
	int status = 0;
	I_T.emplace_back(0);
	for (size_t i = T.size() - 1; i > 0; --i) {
		double T_0 = -T.at(i);
		double T_1 = -T.at(i - 1);
		double h_start = std::fabs(T_0 - T_1);
		//auto d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, h_start, 1e-7, 0.);
		auto d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4imp, h_start, 1e-7, 0.);
		status = gsl_odeiv2_driver_apply(d, &T_0, T_1, Y);
		I_T.emplace_back(Y[0] / -dEdx->get(-T_0));
		gsl_odeiv2_driver_free(d);
	}
	std::reverse(std::begin(I_T), std::end(I_T));
	assert(T.size() == I_T.size());
	return status;
}

int Particle::run_spectrum(const LogAxis& T_) {
	T = T_.get_axis();
	for (auto& T_now : T) {
		I_T.push_back(compute_integral(T_now));
	}
	assert(I_T.size() == T.size());
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

//int Particle::run(const LogAxis& T) {
//	I_T.resize(T.get_size());
//	double Y_n = 0;
//	double r = T.get_ratio();
//	for (double T_n = cgs::TeV; T_n > cgs::GeV; T_n /= 1.01) {
//		//for (auto it = T.get_axis().rbegin(); it != T.get_axis().rend(); ++it) {
//		//	double h_big = *it * (1. - 1. / r);
//		//	double T_n = *it;
//		double h = T_n * (1. - 1. / 1.01);
//		double k_1 = h * f(T_n, Y_n);
//		double k_2 = h * f(T_n - .5 * h, Y_n + .5 * k_1);
//		double k_3 = h * f(T_n - .5 * h, Y_n + .5 * k_2);
//		double k_4 = h * f(T_n - h, Y_n + k_3);
//		Y_n = Y_n + 1. / 6. * (k_1 + 2. * k_2 + 2. * k_3 + k_4);
//		std::cout << T_n / cgs::GeV << " " << Y_n / -dEdx->get(T_n) << "\n";
//
//		//std::cout << std::scientific << *it / cgs::GeV << " " << (*it - h_big) / cgs::GeV << " "
//		//		<< Y_n / -dEdx->get(*it) << "\n";
//	}
//
//	exit(1);
//	return 0;
//}

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
