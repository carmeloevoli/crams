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
	potential = par.potential.get();
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
	return I_Ek.at(i);
}

double Particle::get_modulated(const size_t& i) const {
	return I_Ek_mod.at(i);
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

void Particle::compute_spectrum(const std::vector<double>& Ek) {
	for (auto E_ : Ek) {
		auto I_ = compute_integral(E_);
		I_Ek.push_back(I_);
	}
	assert(I_Ek.size() == Ek.size());
}

void Particle::modulate(const std::vector<double>& Ek) {
	auto E_size = Ek.size();

	double log_E[E_size];
	double log_I[E_size];

	for (int i = 0; i < E_size; ++i) {
		log_E[i] = std::log(Ek.at(i));
		log_I[i] = (I_Ek.at(i) > 0) ? std::log(I_Ek.at(i)) : -100.;
	}

	double phi = std::fabs((double) pid.get_Z() * potential);
	if (pid.get_A())
		phi /= (double) pid.get_A();

	gsl_interp_accel * acc = gsl_interp_accel_alloc();
	gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, E_size);
	gsl_spline_init(spline, log_E, log_I, E_size);

	for (int i = 0; i < E_size; ++i) {
		double factor = 0;
		double T0 = std::min(Ek.at(i) + phi, Ek.back());
		if (pid.get_A())
			factor = (Ek.at(i) * (Ek.at(i) + 2.0 * proton_mass_c2)) / (T0 * (T0 + 2.0 * proton_mass_c2));
		else
			factor = (Ek.at(i) * (Ek.at(i) + 2.0 * electron_mass_c2)) / (T0 * (T0 + 2.0 * electron_mass_c2));
		I_Ek_mod.push_back(factor * std::exp(gsl_spline_eval(spline, std::log(T0), acc)));
	}

	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}

#undef LIMIT
