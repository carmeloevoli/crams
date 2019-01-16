#include <iostream>

#include "axis.h"
#include "cgs.h"
#include "params.h"
#include "particle.h"
#include "output.h"

#include <vector>

int main() {
	Params params;
	params.set_H(4 * cgs::kpc);
	params.print();

	ParticleList particleList;
	particleList.set_abundance(H1, 7e-2);
	particleList.set_abundance(C12, 5e-3);
	particleList.set_abundance(O16, 8e-3);
	particleList.print();

	std::vector<Particle> particles;
	auto list = particleList.get_list();
	for (auto it = list.rbegin(); it != list.rend(); ++it) {
		Particle particle = Particle(it->first, it->second);
		particles.push_back(particle);
	}

	LogAxis T(params.T_min, params.T_max, params.T_size);

	for (auto& particle : particles) {
		std::cout << "running : " << particle.get_pid() << "\n";
		particle.build_grammage(params);
		particle.build_primary_source(params);
		particle.build_inelastic_Xsec(params);
		particle.build_losses(params);
		particle.dump();
		//TODO build secondary_source
		particle.run_spectrum(T);
		particle.clear();
	}

	OutputManager outputManager(T, particles);
	outputManager.dump_spectra(10 * cgs::GeV, 100. * cgs::TeV, 50);
	return 0;
}

