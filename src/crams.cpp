#include "crams.h"

CRAMS::CRAMS(const Params& par_) :
		par(par_) {
	fill_energy(par.E_min.get(), par.E_max.get(), par.E_size.get());
	fill_rigidity(par.E_min.get() * 10, par.E_max.get() * 0.1, par.E_size.get());
	assert(R.size() == T.size());
}

CRAMS::~CRAMS() {
}

void CRAMS::fill_energy(const double& E_min, const double& E_max, const size_t& E_size) {
	for (size_t i = 0; i < E_size; ++i) {
		double ratio = std::pow(E_max / E_min, 1 / (double) (E_size - 1));
		T.push_back(E_min * std::pow(ratio, (double) i));
	}
}

void CRAMS::fill_rigidity(const double& R_min, const double& R_max, const size_t& R_size) {
	for (size_t i = 0; i < R_size; ++i) {
		double ratio = std::pow(R_max / R_min, 1 / (double) (R_size - 1));
		R.push_back(R_min * std::pow(ratio, (double) i));
	}
}

void CRAMS::fill_particles(const NucleiList& nucleilist) {
	for (auto nucleus : nucleilist.list) { // TODO pass by reference this https://stackoverflow.com/questions/12851860/accessing-a-map-by-returning-a-pointer-to-it-c
		Particle particle = Particle(nucleus.first, nucleus.second.first, nucleus.second.second, par);
		particles.push_back(particle);
	}
	assert(particles.size() == nucleilist.get_size());
}

void CRAMS::run() {
	for (auto &p : particles) {
		std::cout << "# " << p.get_pid() << "\n";
		p.compute_spectrum(T);
		p.modulate(T, R);
	}
}

void CRAMS::dump() {
	std::cout << "# T\t";
	for (auto particle : particles)
		std::cout << particle.get_pid() << "\t";
	std::cout << "\n";
	for (size_t i = 0; i < T.size(); ++i) {
		std::cout << std::scientific << T.at(i) / GeV << "\t";
		for (auto particle : particles)
			//std::cout << p.get_modulated(i) / (1. / GeV / m2 / s / sr) << "\t";
			std::cout << particle.get_I(i) / (1. / GeV / m2 / sec / sr) << "\t";
		std::cout << "\n";
	}
}
