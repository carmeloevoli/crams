#include "output.h"

#include <plog/Log.h>

#include <iomanip>
#include <iostream>

#include "utilities.h"

namespace CRAMS {

OutputManager::OutputManager(const Particles& particles, const Input& input)
    : m_particles(particles),  // TODO is the best way to pass particle?
      m_phi(input.modulationPotential),
      m_id(input.id) {
  m_R = Utilities::LogAxis(input.ROutputMin, input.ROutputMax, input.ROutputSize);
}

OutputManager::~OutputManager() { LOGD << "deleted OutputManager"; }

// ptr_Particle OutputManager::find_ptr(const PID& pid) {
//   auto it = find(_particles.begin(), _particles.end(), Particle(pid, 0));
//   bool isPresent = !(it == _particles.end());
//   return {isPresent, it};
// }

double OutputManager::getFluxChargeGroup(const int Z, const double& R) const {
  double value = 0.;
  for (const auto& particle : m_particles) {
    if (particle.isChargeZ(Z)) {
      value += particle.I_R_TOA(R, m_phi);
    }
  }
  return value;
}

void OutputManager::dumpSpectra() const {
  const std::string spectraFilename = "output/crams_spectra_" + std::to_string(m_id) + ".txt";
  LOGW << "writing spectra to " << spectraFilename;
  std::ofstream outfile(spectraFilename.c_str());
  if (outfile.is_open()) {
    outfile << std::scientific;
    const double units = 1. / (CGS::GeV * pow2(CGS::meter) * CGS::sec);
    for (const auto& R_i : m_R) {
      outfile << R_i / CGS::GeV << "\t";
      for (int iZ = 1; iZ <= 28; ++iZ) outfile << getFluxChargeGroup(iZ, R_i) / units << "\t";
      outfile << "\n";
    }
    outfile.close();
  } else {
    throw std::runtime_error("file not open for writing");
  }
}

}  // namespace CRAMS