#include <iostream>

#include "cgs.h"
#include "params.h"
#include "particle.h"
#include "output.h"
#include "utilities.h"

#include <vector>

int main() {
	//TODO read params from ini file
	//TODO print the commit number

	Params params;
	params.set_H(4 * cgs::kpc);
	params.print();

	ParticleList particleList;
	particleList.set_abundance(H1, 7e-2);
	particleList.set_abundance(C12, 5e-3);
	particleList.set_abundance(O16, 8e-3);
	particleList.set_abundance(Ne20, 1e-4);
	particleList.set_abundance(Mg24, 1e-4);
	particleList.set_abundance(Si28, 1e-4);
	particleList.set_abundance(Fe56, 1e-4);
	particleList.print();

	std::vector<Particle> particles;
	auto list = particleList.get_list();
	for (auto it = list.rbegin(); it != list.rend(); ++it) {
		Particle particle = Particle(it->first, it->second);
		particles.push_back(particle);
	}

	auto T = LogAxis(params.T_min, params.T_max, params.T_size);

	for (auto& particle : particles) {
		std::cout << "running : " << particle.get_pid() << "\n";
		particle.build_grammage(params);
		particle.build_primary_source(params);
		particle.build_secondary_source(particles);
		particle.build_inelastic_Xsec();
		particle.build_losses(params);
		//particle.dump();
		particle.setDone() = particle.run(T);
		particle.clear();
	}

	OutputManager outputManager(particles, params.modulation_potential);
	outputManager.dump_spectra(10 * cgs::GeV, 10. * cgs::TeV, 50);
	outputManager.dump_heavy_spectra(10 * cgs::GeV, 10. * cgs::TeV, 50);
	outputManager.dump_ratio(10 * cgs::GeV, 10. * cgs::TeV, 50);
	return 0;
}

