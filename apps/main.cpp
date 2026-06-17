#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "crams.h"

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

    std::unique_ptr<CRAMS::InelasticXsec> inelasticXsecs;
    switch (input.inelasticModel()) {
      case CRAMS::InelasticModel::Glauber:
        inelasticXsecs = std::make_unique<CRAMS::InXsecGlauber>();
        break;
      case CRAMS::InelasticModel::Tripathi99:
        inelasticXsecs = std::make_unique<CRAMS::InXsecTripathi99>();
        break;
      case CRAMS::InelasticModel::Crosec:
        inelasticXsecs = std::make_unique<CRAMS::InXsecCrosec>();
        break;
    }

    std::unique_ptr<CRAMS::NucFragXsec> nucfragXsecs;
    switch (input.fragmentationModel()) {
      case CRAMS::FragmentationModel::Fluka4Dragon:
        nucfragXsecs = std::make_unique<CRAMS::NucFragFluka4Dragon>();
        break;
      case CRAMS::FragmentationModel::UsineGalprop17Opt12:
        nucfragXsecs = std::make_unique<CRAMS::NucFragUsineGalprop17Opt12>();
        break;
      case CRAMS::FragmentationModel::UsineGalprop17Opt22:
        nucfragXsecs = std::make_unique<CRAMS::NucFragUsineGalprop17Opt22>();
        break;
      case CRAMS::FragmentationModel::UsineWebber03Coste12:
        nucfragXsecs = std::make_unique<CRAMS::NucFragUsineWebber03Coste12>();
        break;
      case CRAMS::FragmentationModel::Evoli2026W93:
        nucfragXsecs = std::make_unique<CRAMS::NucFragEvoli2026W93>();
        break;
      case CRAMS::FragmentationModel::Evoli2026St99:
        nucfragXsecs = std::make_unique<CRAMS::NucFragEvoli2026St99>();
        break;
    }

    for (auto& particle : particles) {
      LOGI << "running : " << particle.getPid();
      particle.buildVectors(input);
      particle.buildGrammage(input);
      particle.buildLosses(input);
      particle.buildPrimarySource(input);
      particle.buildInelasticXsecs(*inelasticXsecs);
      if (!particle.getPid().isTertiary()) particle.buildSecondarySource(input, particles, *nucfragXsecs);
      if (particle.getPid() == CRAMS::H1_ter) particle.buildTertiarySource(particles);
      //  if (input.X_s() > 0.) particle.buildGrammageAtSource(input, particles, *nucfragXsecs);
      if (!quiet) particle.dump();
      particle.computeIntensity(input);
      particle.reset();
    }

    CRAMS::OutputManager outputManager(particles, input);
    outputManager.dumpSpectraRigidity();
    outputManager.dumpIsotopes();  // cheap (2 columns); needed by the MCMC Be10/Be9 posterior
    if (!quiet) {
      outputManager.dumpSpectraEkn();
    }
  } catch (const std::exception& e) {
    LOGE << "exception caught with message: " << e.what();
  }
  return 0;
}
