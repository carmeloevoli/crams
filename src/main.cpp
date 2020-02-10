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
#ifdef DEBUG
	std::cout << "version: " << get_version() << "\n";
	std::cout << "git version: " << git_sha1() << "\n";
	std::cout << "has local changes: " << std::boolalpha << git_has_local_changes()
	<< std::noboolalpha << "\n";
	std::cout << "was built on: " << __DATE__ << " " __TIME__ << "\n";
#endif
}

int main(int argc, char * argv[]) {
	if (argc == 2) {
		log_startup_information();

		Params params;
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
#ifdef DEBUG
			std::cout << "running : " << particle.get_pid() << "\n";
#endif
			particle.build_grammage(params);
			particle.build_snr_source(params);
			particle.build_inelastic_Xsec(params);
			particle.build_secondary_source(particles, params);
			if (particle.get_pid() == H1_ter)
				particle.build_tertiary_source(particles);
//			particle.build_grammage_at_source(particles, params);
			particle.build_losses(params);
//			particle.dump();
			particle.setDone() = particle.run(T);
			particle.clear();
		}

		OutputManager outputManager(particles, params.modulation_potential, params.id);
	        outputManager.dump_spectra(2 * cgs::GeV, 5. * cgs::TeV, 100);
		outputManager.dump_ratios(2 * cgs::GeV, 5. * cgs::TeV, 100);
//      outputManager.dump_heavy_spectra(10 * cgs::GeV, 10. * cgs::TeV, 50);

		std::ofstream fchi2("chi2_results.txt", std::ofstream::out);

		Chi2_H chi2_H(particles, params.modulation_potential);
		fchi2 << chi2_H.compute_chi2(10. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_He chi2_He(particles, params.modulation_potential);
		fchi2 << chi2_He.compute_chi2(10. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_C chi2_C(particles, params.modulation_potential);
		fchi2 << chi2_C.compute_chi2(10. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_N chi2_N(particles, params.modulation_potential);
		fchi2 << chi2_N.compute_chi2(10. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_O chi2_O(particles, params.modulation_potential);
		fchi2 << chi2_O.compute_chi2(10. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_HeO chi2_HeO(particles, params.modulation_potential);
		fchi2 << chi2_HeO.compute_chi2(2. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_BeB chi2_BeB(particles, params.modulation_potential);
		fchi2 << chi2_BeB.compute_chi2(2. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_BeC chi2_BeC(particles, params.modulation_potential);
		fchi2 << chi2_BeC.compute_chi2(2. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_BeO chi2_BeO(particles, params.modulation_potential);
		fchi2 << chi2_BeO.compute_chi2(2. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_BC chi2_BC(particles, params.modulation_potential);
		fchi2 << chi2_BC.compute_chi2(2. * cgs::GeV, 500 * cgs::GeV) << "\n";

		Chi2_BO chi2_BO(particles, params.modulation_potential);
		fchi2 << chi2_BO.compute_chi2(2. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_CO chi2_CO(particles, params.modulation_potential);
		fchi2 << chi2_CO.compute_chi2(2. * cgs::GeV, 500. * cgs::GeV) << "\n";

		Chi2_BeB_statsonly chi2_BeB_statsonly(particles, params.modulation_potential);
		fchi2 << chi2_BeB_statsonly.compute_chi2(2. * cgs::GeV, 500. * cgs::GeV) << "\n";

		fchi2.close();

	} else {
		std::cout << "Usage: ./CRAMS params.ini\n";
	}
	return 0;
}

