#include <iostream>
#include <vector>

#include "crams/core/cgs.h"
#include "crams/core/input.h"
#include "crams/core/output.h"
#include "crams/particle.h"
#include "crams/utils/logging.h"
#include "crams/utils/utilities.h"

int main(int argc, char* argv[]) {
  log_startup_information();
  try {
    CRAMS::Input input;
    CRAMS::ParticleList particleList;

    if (argc == 2) {
      input.setSimname(argv[1]);
      input.readParamsFromFile(argv[1]);
      input.print();
      particleList.readParamsFromFile(argv[1]);
      particleList.print();
    } else if (argc == 1) {
      LOGI << "no input file provided, using default parameters";
    } else {
      throw std::runtime_error("too many arguments provided, expected './crams params.ini'");
    }

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
      // particle.buildInelasticXsecs(input);
      //  particle.buildSecondarySource(input, particles);
      //  if (particle.getPid() == CRAMS::H1_ter) particle.buildTertiarySource(particles);
      //  if (input.X_s() > 0.) particle.buildGrammageAtSource(input, particles);
      particle.dump();
      particle.computeIntensity(input);
      particle.reset();
    }

    CRAMS::OutputManager outputManager(particles, input);
    outputManager.dumpSpectraRigidity();
    outputManager.dumpSpectraEkn();
    outputManager.dumpIsotopes();
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}
