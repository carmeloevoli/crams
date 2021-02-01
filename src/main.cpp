#include <iostream>
#include <vector>

#include "cgs.h"
// #include "chi2.h"
#include "git_revision.h"
#include "input.h"
#include "logging.h"
#include "output.h"
#include "particle.h"
#include "utilities.h"

int main(int argc, char* argv[]) {
  log_startup_information();
  try {
    if (argc == 2) {
      CRAMS::Input input;
      input.readParamsFromFile(argv[1]);
      input.print();

      CRAMS::ParticleList particleList;
      particleList.readParamsFromFile(argv[1]);
      particleList.print();

      CRAMS::Particles particles;
      auto list = particleList.getList();
      particles.reserve(list.size());
      for (auto it = list.rbegin(); it != list.rend(); ++it) {
        auto pid = it->first;
        auto nucleusParams = it->second;
        particles.emplace_back(pid, nucleusParams);
      }

      for (auto& particle : particles) {
        LOGI << "running : " << particle.getPid();
        particle.buildVectors(input);
        particle.buildGrammage(input);
        particle.buildLosses(input);
        particle.buildPrimarySource(input);
        particle.buildInelasticXsecs(input);
        particle.buildSecondarySource(input, particles);
        //         // // 			if (particle.get_pid() == H1_ter)
        //         // // 				particle.build_tertiary_source(particles);
        //         if (input.X_s > 0.) particle.buildGrammageAtSource(input, particles);
        particle.computeIntensity();
        particle.dump();
        particle.reset();
      }

      CRAMS::OutputManager outputManager(particles, input);
      outputManager.dumpSpectra();

      // 		std::ofstream fchi2("chi2_results.txt", std::ofstream::out);

      // 		double T_min = 3. * cgs::GeV;
      // 		double T_max = 500. * cgs::GeV;

      // 		Chi2_H chi2_H(particles, params.modulation_potential);
      // 		fchi2 << chi2_H.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_He chi2_He(particles, params.modulation_potential);
      // 		fchi2 << chi2_He.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_B chi2_B(particles, params.modulation_potential);
      // 		fchi2 << chi2_B.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_C chi2_C(particles, params.modulation_potential);
      // 		fchi2 << chi2_C.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_N chi2_N(particles, params.modulation_potential);
      // 		fchi2 << chi2_N.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_O chi2_O(particles, params.modulation_potential);
      // 		fchi2 << chi2_O.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_HeO chi2_HeO(particles, params.modulation_potential);
      // 		fchi2 << chi2_HeO.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_BeB chi2_BeB(particles, params.modulation_potential);
      // 		fchi2 << chi2_BeB.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_BeC chi2_BeC(particles, params.modulation_potential);
      // 		fchi2 << chi2_BeC.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_BeO chi2_BeO(particles, params.modulation_potential);
      // 		fchi2 << chi2_BeO.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_BC chi2_BC(particles, params.modulation_potential);
      // 		fchi2 << chi2_BC.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_BO chi2_BO(particles, params.modulation_potential);
      // 		fchi2 << chi2_BO.compute_chi2(T_min, T_max) << "\n";

      // 		Chi2_CO chi2_CO(particles, params.modulation_potential);
      // 		fchi2 << chi2_CO.compute_chi2(T_min, T_max) << "\n";

      // 		//Chi2_BeB_statsonly chi2_BeB_statsonly(particles, params.modulation_potential);
      // 		//fchi2 << chi2_BeB_statsonly.compute_chi2(2. * cgs::GeV, 500. * cgs::GeV) << "\n";

      // 		fchi2.close();
    } else {
      throw std::runtime_error("please provide an input file as './crams params.ini'");
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}
