
#include "crams/runner.h"

#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

CRAMS::Runner::Runner(InelasticModel inelasticModel, FragmentationModel fragmentationModel, ParticleList injection)
    : inelasticModel{inelasticModel}, fragmentationModel{fragmentationModel}, injection{injection} {
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

void CRAMS::Runner::setInjectionParams(std::vector<double> abundances, std::vector<double> slopes) {
  auto Zmin = injection.lightest().getZ();
  auto Zmax = injection.heaviest().getZ();

  for (std::size_t abIdx = 0; abIdx < abundances.size(); ++abIdx) {
    injection.setAbundanceChargeGroup(Zmin + abIdx, abundances[abIdx]);
  }

  auto allSlopesSpecified = slopes.size() >= (Zmax - Zmin + 1);
  auto endSlopeIdx = allSlopesSpecified ? slopes.size() : slopes.size() - 1;
  for (std::size_t slopeIdx = 0; slopeIdx < endSlopeIdx; ++slopeIdx) {
    injection.setSlopeChargeGroup(Zmin + slopeIdx, slopes[slopeIdx]);
  }
  if (!allSlopesSpecified) {
    injection.setSlopeNuclei(Zmin + endSlopeIdx, slopes.back());
  }
}

CRAMS::RigiditySpectra CRAMS::Runner::compute(Input input, bool dumpToFile, bool verbose, bool ignoreInputInitParams) {
  if (!ignoreInputInitParams &&
      ((input.inelasticModel() != inelasticModel) || (input.fragmentationModel() != fragmentationModel))) {
    throw std::runtime_error("compute method called on input with mismatching inelastic and/or fragmentation model");
  }

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
    if (verbose && dumpToFile) particle.dump();
    particle.computeIntensity(input);
    particle.reset();
  }

  OutputManager outputManager(result, input);
  if (dumpToFile) {
    outputManager.dumpSpectraRigidity();
    outputManager.dumpIsotopes();  // cheap (2 columns); needed by the MCMC Be10/Be9 posterior
    if (verbose) {
      outputManager.dumpSpectraEkn();
    }
  }

  return outputManager.rigiditySpectra();
}
