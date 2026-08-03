
#include "crams/runner.h"

#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

CRAMS::Runner::Runner(InelasticModel im, FragmentationModel fm) {
  inelasticModel = im;
  fragmentationModel = fm;

  switch (inelasticModel) {
    case CRAMS::InelasticModel::Tripathi99:
      inelasticXsecs = std::make_unique<CRAMS::InXsecTripathi99>();
      break;
    case CRAMS::InelasticModel::Glauber:
      inelasticXsecs = std::make_unique<CRAMS::InXsecGlauber>();
      break;
    case CRAMS::InelasticModel::Crosec:
      inelasticXsecs = std::make_unique<CRAMS::InXsecCrosec>();
      break;
  }

  switch (fragmentationModel) {
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
    case CRAMS::FragmentationModel::Evoli2019:
      nucfragXsecs = std::make_unique<CRAMS::NucFragEvoli2019>();
      break;
    case CRAMS::FragmentationModel::Evoli2026W93:
      nucfragXsecs = std::make_unique<CRAMS::NucFragEvoli2026W93>();
      break;
    case CRAMS::FragmentationModel::Evoli2026St99:
      nucfragXsecs = std::make_unique<CRAMS::NucFragEvoli2026St99>();
      break;
  }
}

CRAMS::Particles CRAMS::Runner::compute(ParticleList injection, Input input, bool dumpToFile, bool verbose) {
  Particles result;
  auto list = injection.getList();
  result.reserve(list.size());
  for (auto it = list.rbegin(); it != list.rend(); ++it) {
    auto pid = it->first;
    auto nucleusParams = it->second;
    result.emplace_back(pid, nucleusParams);
  }

  for (auto& particle : result) {
    if (verbose) {
      LOGI << "running : " << particle.getPid();
    }
    particle.buildVectors(input);
    particle.buildGrammage(input);
    particle.buildLosses(input);
    particle.buildPrimarySource(input);
    particle.buildInelasticXsecs(*inelasticXsecs);
    if (!particle.getPid().isTertiary()) particle.buildSecondarySource(input, result, *nucfragXsecs);
    if (particle.getPid() == CRAMS::H1_ter) particle.buildTertiarySource(result);
    //  if (input.X_s() > 0.) particle.buildGrammageAtSource(input, result, *nucfragXsecs);
    if (verbose) particle.dump();
    particle.computeIntensity(input);
    particle.reset();
  }

  if (dumpToFile) {
    OutputManager outputManager(result, input);
    outputManager.dumpSpectraRigidity();
    outputManager.dumpIsotopes();  // cheap (2 columns); needed by the MCMC Be10/Be9 posterior
    if (verbose) {
      outputManager.dumpSpectraEkn();
    }
  }

  return result;
}
