#include <iostream>
#include <vector>

#include "cgs.h"
#include "chi2.h"
#include "params.h"
#include "particle.h"
#include "output.h"
#include "utilities.h"

#include "git_revision.h"

void log_startup_information() {
	std::cout << "version: " << get_version() << "\n";
	std::cout << "git version: " << git_sha1() << "\n";
	std::cout << "has local changes: " << std::boolalpha << git_has_local_changes()
			<< std::noboolalpha << "\n";
	std::cout << "was built on: " << __DATE__ << " " __TIME__ << "\n";
}

int main(int argc, char * argv[]) {
	if (argc == 2) {
		log_startup_information();

		Params params;
		params.set_H(4 * cgs::kpc);
		params.set_D0(1.8e28 * cgs::cm2 / cgs::sec);
		params.set_from_file(argv[1]);
		params.print();

		ParticleList particleList;
		particleList.set_from_file(argv[1]);
		particleList.print();

		Particles particles;
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
		//outputManager.dump_heavy_spectra(10 * cgs::GeV, 10. * cgs::TeV, 50);
		outputManager.dump_ratio(10 * cgs::GeV, 10. * cgs::TeV, 50);

		std::ofstream fchi2("chi2_results.txt", std::ofstream::out);

		Chi2_C chi2_C(particles, params.modulation_potential);
		fchi2 << chi2_C.compute_chi2(20. * cgs::GeV, 1e3 * cgs::GeV) << "\n";

		Chi2_N chi2_N(particles, params.modulation_potential);
		fchi2 << chi2_N.compute_chi2(20. * cgs::GeV, 300. * cgs::GeV) << "\n";

		Chi2_O chi2_O(particles, params.modulation_potential);
		fchi2 << chi2_O.compute_chi2(20. * cgs::GeV, 1e3 * cgs::GeV) << "\n";

		Chi2_BC chi2_BC(particles, params.modulation_potential);
		fchi2 << chi2_BC.compute_chi2(20. * cgs::GeV, 300. * cgs::GeV) << "\n";

		Chi2_CO chi2_CO(particles, params.modulation_potential);
		fchi2 << chi2_CO.compute_chi2(20. * cgs::GeV, 300. * cgs::GeV) << "\n";

		fchi2.close();

	} else {
		std::cout << "Usage: ./CRAMS params.ini\n";
	}
	return 0;
}

