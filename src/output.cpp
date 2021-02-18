#include "output.h"

#include <plog/Log.h>

#include <iomanip>
#include <iostream>

#include "utilities.h"

namespace CRAMS {

OutputManager::OutputManager(const Particles& particles, const Input& input)
    : m_particles(particles),  // TODO is the best way to pass particle?
      m_phi(input.modulationPotential),
      m_id(input.id),
      m_simname(input.simname) {
  m_R = Utilities::LogAxis(input.ROutputMin, input.ROutputMax, input.ROutputSize);
}

OutputManager::~OutputManager() { LOGD << "deleted OutputManager"; }

// ptr_Particle OutputManager::find_ptr(const PID& pid) {
//   auto it = find(_particles.begin(), _particles.end(), Particle(pid, 0));
//   bool isPresent = !(it == _particles.end());
//   return {isPresent, it};
// }

double OutputManager::I_R_TOA(const Particle& particle, const double& R, const double& modulationPotential) const {
  // see arXiv:1511.08790
  const auto pid = particle.getPid();
  double value = 0;
  {
    constexpr double mpSquared = pow2(CGS::protonMassC2);
    const double ZOverASquared = pow2(pid.getZoverA());
    const double ESquared = pow2(R) * ZOverASquared + mpSquared;
    const double T = std::sqrt(ESquared) - CGS::protonMassC2;
    const double Phi = pid.getZoverA() * modulationPotential;
    const double T_ISM = T + Phi;
    double dTdR = R * ZOverASquared;
    dTdR /= std::sqrt(ZOverASquared * pow2(R) + mpSquared);
    double factor = T * (T + 2. * CGS::protonMassC2);
    factor /= (T + Phi) * (T + Phi + 2. * CGS::protonMassC2);
    value = factor * particle.I_T_interpol(T_ISM) * dTdR;
  }
  return value;
}

double OutputManager::getFluxChargeGroup(const int Z, const double& R) const {
  double value = 0.;
  for (const auto& particle : m_particles) {
    if (particle.isChargeZ(Z)) {
      value += I_R_TOA(particle, R, m_phi);
    }
  }
  return value;
}

void OutputManager::dumpSpectra() const {
  const std::string spectraFilename = "output/" + m_simname + "_spectra_R_" + std::to_string(m_id) + ".txt";
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

double OutputManager::getFluxChargeGroupEkn(const int Z, const double& T) const {
  double value = 0.;
  for (const auto& particle : m_particles) {
    if (particle.isChargeZ(Z)) {
      const double Phi = m_phi * particle.getPid().getZoverA();
      double factor = T * (T + 2. * CGS::protonMassC2);
      factor /= (T + Phi) * (T + Phi + 2. * CGS::protonMassC2);
      value += factor * particle.I_T_interpol(T + Phi);
    }
  }
  return value;
}

void OutputManager::dumpSpectraEkn() const {
  const std::string spectraFilename = "output/" + m_simname + "_spectra_Ekn_" + std::to_string(m_id) + ".txt";
  LOGW << "writing spectra to " << spectraFilename;
  std::ofstream outfile(spectraFilename.c_str());
  auto T = Utilities::LogAxis(0.3 * CGS::GeV, 1. * CGS::TeV, 4 * 32);
  const double units = 1. / (CGS::GeV * pow2(CGS::meter) * CGS::sec);
  if (outfile.is_open()) {
    outfile << std::scientific;
    for (const auto& T_i : T) {
      outfile << T_i / CGS::GeV << "\t";
      for (int iZ = 1; iZ <= 28; ++iZ) outfile << getFluxChargeGroupEkn(iZ, T_i) / units << "\t";
      outfile << "\n";
    }
    outfile.close();
  } else {
    throw std::runtime_error("file not open for writing");
  }
}

}  // namespace CRAMS