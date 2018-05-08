#include "crams.h"

CRAMS::CRAMS(const Params& par_) :
		par(par_) {
	fill_energy(par.E_min.get(), par.E_max.get(), par.E_size.get());
	add_nucleus(PID(6, 12), 0.01, 4.2);
	add_nucleus(PID(8, 16), 0.01, 4.2);
	fill_particles();
}

CRAMS::~CRAMS() {
}

void CRAMS::fill_energy(const double& E_min, const double& E_max, const size_t& E_size) {
	for (size_t i = 0; i < E_size; ++i) {
		double ratio = std::pow(E_max / E_min, 1 / (double) (E_size - 1));
		E.push_back(E_min * std::pow(ratio, (double) i));
	}
}

void CRAMS::add_nucleus(const PID& pid, const double& efficiency, const double& gamma) {
	std::pair<double, double> p = std::make_pair(efficiency, gamma);
	nucleilist[pid] = p;
}

void CRAMS::fill_particles() {
	for (auto nucleus : nucleilist) {
		Particle particle = Particle(nucleus.first, nucleus.second.first, nucleus.second.second, par);
		particles.push_back(particle);
	}
}

void CRAMS::test() {
	std::cout << particles.size() << "\n";
	particles[0].dump_inputs();
}

void CRAMS::run() {
	for (auto &p : particles) {
		std::cout << "# " << p.get_pid() << "\n";
		p.compute_spectrum(E);
		p.modulate(E);
	}
}

void CRAMS::dump() {
	std::cout << "# E_k\t";
	for (auto p : particles)
		std::cout << p.get_pid() << "\t";
	std::cout << "\n";
	for (size_t i = 0; i < E.size(); ++i) {
		std::cout << std::scientific << E.at(i) / GeV << "\t";
		for (auto p : particles)
			std::cout << p.get_I(i) << "\t" << p.get_modulated(i) << "\t";
		std::cout << "\n";
	}
}
