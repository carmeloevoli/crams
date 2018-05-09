#include "particle.h"

#define LIMIT 1000

Particle::Particle() {
}

Particle::Particle(PID pid_, double efficiency_, double slope_, Params par) :
		pid(pid_) {
	X = Grammage(pid, par);
	Q = PrimarySource(pid, par.mu.get(), efficiency_, slope_);
	sigma_in = InelasticXsec(pid);
	b = Losses(pid, par);
	modulation_potential = par.potential.get();
}

Particle::~Particle() {
}

void Particle::dump_inputs() const {
	std::cout << std::scientific;
	for (double E = 0.1 * GeV; E < 10 * PeV; E *= 1.1) {
		std::cout << E / GeV << "\t";
		std::cout << X.get(E) / (gram / cm2) << "\t";
		std::cout << X.diffusion_escape_time(E) / year << "\t";
		std::cout << X.advection_escape_time() / year << "\t";
		std::cout << std::pow(E, 2.2) * Q.get(E) << "\t";
		std::cout << -(E / b.dE_dt_adiabatic(E)) / year << "\t";
		std::cout << -(E / b.dE_dt_ionization(E)) / year << "\t";
		std::cout << sigma_in.get(E) / mbarn << "\t";
		std::cout << "\n";
	}
}

PID Particle::get_pid() const {
	return pid;
}

double Particle::get_I(const size_t& i) const {
	return I_T.at(i);
}

double Particle::get_modulated(const size_t& i) const {
	return I_R_TOA.at(i);
}

double compute_integral_qags(gsl_integration_workspace * w, gsl_function * F, double x_lo, double x_hi) {
	double result, error;
	gsl_integration_qags(F, x_lo, x_hi, 0, 1e-4, LIMIT, w, &result, &error);
	return result;
}

double Particle::internal_integrand(const double& E_second) {
	return Lambda_1(E_second) / Lambda_2(E_second);
}

double gslParticleClassExpWrapper(double x, void * pp) {
	double E_second = std::exp(x);
	Particle * particle = (Particle *) pp;
	return E_second * particle->internal_integrand(E_second);
}

double Particle::ExpIntegral(const double& E, const double& E_prime) {
	double result = 0;

	gsl_function F;
	F.params = this;
	F.function = &gslParticleClassExpWrapper;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(LIMIT);

	result = compute_integral_qags(w, &F, std::log(E), std::log(E_prime));

	gsl_integration_workspace_free(w);

	return result;
}

double Particle::external_integrand(const double& E_prime, const double& E) {
	return Q.get(E_prime) / Lambda_2(E_prime) * std::exp(-ExpIntegral(E, E_prime));
}

double Particle::Lambda_1(const double& E) {
	double value = 1. / X.get(E) + sigma_in.get(E) / mean_ism_mass - b.get_derivative(E);
	return value;
}

double Particle::Lambda_2(const double& E) {
	double value = b.dE_dx(E);
	return -value;
}

struct gsl_f_pars {
	double E;
	Particle * pt_Particle;
};

double gslParticleClassWrapper(double x, void * pp) {
	gsl_f_pars *p = (gsl_f_pars *) pp;
	double E = p->E;
	double E_prime = std::exp(x);
	return E_prime * p->pt_Particle->external_integrand(E_prime, E);
}

double Particle::compute_integral(const double& E) {
	double result = 0;

	gsl_f_pars pars = { E, this };
	gsl_function F;
	F.params = &pars;
	F.function = &gslParticleClassWrapper;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(LIMIT);

	result = compute_integral_qags(w, &F, std::log(E), std::log(E) + 2);

	gsl_integration_workspace_free(w);

	return result;
}

void Particle::compute_spectrum(const std::vector<double>& T) {
	for (auto T_ : T) {
		auto I_ = compute_integral(T_);
		I_T.push_back(I_);
	}
	assert(I_T.size() == T.size());
}

void Particle::modulate(const std::vector<double>& T, const std::vector<double>& R) {
	auto T_size = T.size();
	double log_T[T_size];
	double log_I[T_size];

	for (int i = 0; i < T_size; ++i) {
		log_T[i] = std::log(T.at(i));
		log_I[i] = (I_T.at(i) > 0) ? std::log(I_T.at(i)) : -100.;
	}

	double Phi = pid.get_Z_over_A() * modulation_potential;

	gsl_interp_accel * acc = gsl_interp_accel_alloc();
	gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, T_size);
	gsl_spline_init(spline, log_T, log_I, T_size);

	for (int i = 0; i < R.size(); ++i) {
		double T_ = std::sqrt(pow2(pid.get_Z_over_A() * R.at(i)) + pow2(proton_mass_c2)) - proton_mass_c2;
		const double T_phi = std::min(T_ + Phi, T.back());
		double factor = pid.get_Z_over_A();
		factor *= std::sqrt(T_ * (T_ + 2. * proton_mass_c2)) / (T_ + proton_mass_c2);
		factor *= T_ * (T_ + 2.0 * proton_mass_c2) / T_phi / (T_phi + 2.0 * proton_mass_c2);
		I_R_TOA.push_back(factor * std::exp(gsl_spline_eval(spline, std::log(T_phi), acc)));
	}

	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}

#undef LIMIT
