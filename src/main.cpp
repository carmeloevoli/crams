#include <iostream>
#include <vector>

#include "cgs.h"
#include "chi2.h"
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
      input.setSimname(argv[1]);
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
        if (particle.getPid() == CRAMS::H1_ter) particle.buildTertiarySource(particles);
        if (particle.getPid() == CRAMS::pbar) particle.buildAntiprotonSource(particles);
        // if (input.X_s > 0.) particle.buildGrammageAtSource(input, particles);
        particle.dump();
        particle.computeIntensity();
        particle.reset();
      }

      CRAMS::OutputManager outputManager(particles, input);
      outputManager.dumpSpectra();
      outputManager.dumpSpectraEkn();

      {
        std::ofstream fchi2("chi2_results.txt", std::ofstream::out);
        const std::pair<double, double> R_range = {20. * CRAMS::CGS::GeV, 400. * CRAMS::CGS::GeV};
        // std::vector<std::unique_ptr<CRAMS::Chi2>> chi2s;

        // chi2s.emplace_back(new CRAMS::Chi2IH("H", particles, input.modulationPotential));
        // chi2s.emplace_back(new CRAMS::Chi2IHe("He", particles, input.modulationPotential));

        // for (const auto& c : chi2s) {
        //   fchi2 << c->getName() << ", " << c->computeChi2(R_range.first, R_range.second) << "\n";
        // }
        fchi2.close();
      }
    } else {
      throw std::runtime_error("please provide an input file as './crams params.ini'");
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}
