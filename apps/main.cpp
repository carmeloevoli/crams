#include <iostream>
#include <vector>

#include "crams/core/cgs.h"
#include "crams/core/input.h"
#include "crams/core/output.h"
#include "crams/inelastic.h"
#include "crams/particle.h"
#include "crams/particlelist.h"
#include "crams/utils/logging.h"
#include "crams/utils/utilities.h"

int main(int argc, char* argv[]) {
  bool quiet = false;
  std::string inifile;

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "-q" || arg == "--quiet") {
      quiet = true;
    } else if (inifile.empty()) {
      inifile = arg;
    } else {
      throw std::runtime_error("unexpected argument: '" + arg + "'");
    }
  }

  log_startup_information(quiet);

  try {
    CRAMS::Input input;
    CRAMS::ParticleList particleList;

    if (!inifile.empty()) {
      input.setSimname(inifile);
      input.readParamsFromFile(inifile);
      if (!quiet) input.print();
      particleList.readParamsFromFile(inifile);
      if (!quiet) particleList.print();
    } else {
      LOGI << "no input file provided, using default parameters";
    }

    CRAMS::Particles particles;
    auto list = particleList.getList();
    particles.reserve(list.size());
    for (auto it = list.rbegin(); it != list.rend(); ++it) {
      auto pid = it->first;
      auto nucleusParams = it->second;
      particles.emplace_back(pid, nucleusParams);
    }

    CRAMS::InXsecTripathi99 inelasticXsecs;
    for (auto& particle : particles) {
      LOGI << "running : " << particle.getPid();
      particle.buildVectors(input);
      particle.buildGrammage(input);
      particle.buildLosses(input);
      particle.buildPrimarySource(input);
      particle.buildInelasticXsecs(inelasticXsecs);
      //  particle.buildSecondarySource(input, particles);
      //  if (particle.getPid() == CRAMS::H1_ter) particle.buildTertiarySource(particles);
      //  if (input.X_s() > 0.) particle.buildGrammageAtSource(input, particles);
      particle.dump();
      particle.computeIntensity(input);
      particle.reset();
    }

    CRAMS::OutputManager outputManager(particles, input);
    outputManager.dumpSpectraRigidity();
    if (!quiet) {
      outputManager.dumpSpectraEkn();
      outputManager.dumpIsotopes();
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}
