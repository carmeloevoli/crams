#include "crams/core/output.h"

#include <plog/Log.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "crams/utils/utilities.h"

namespace CRAMS {

using Utilities::pow2;

namespace {

const std::array<const char*, 29> kElementSymbols = {{
    "",  "H", "He", "Li", "Be", "B",  "C",  "N",  "O", "F",  "Ne", "Na", "Mg", "Al", "Si",
    "P", "S", "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
}};

std::string elementSymbol(int Z) {
  if (Z > 0 && static_cast<size_t>(Z) < kElementSymbols.size()) return kElementSymbols[Z];
  return "Z" + std::to_string(Z);
}

void writeFluxUnit(std::ostream& out) { out << "# Flux unit -> 1 / (GeV m2 s sr)\n"; }

void writeChargeGroupColumns(std::ostream& out, const std::string& axisLabel) {
  out << "# Columns\n";
  out << "# " << axisLabel << " -> 1\n";
  out << "# Flux columns are charge groups summed over isotopes.\n";

  int column = 2;
  for (int Z = 1; Z <= 28; ++Z, ++column) {
    out << "# " << elementSymbol(Z) << " -> " << column << "\n";
  }
  writeFluxUnit(out);
}

void writeIsotopeColumns(std::ostream& out) {
  out << "# Columns\n";
  out << "# R [GV] -> 1\n";
  out << "# Be7 -> 2\n";
  out << "# Be9 -> 3\n";
  out << "# Be10 -> 4\n";
  writeFluxUnit(out);
}

}  // namespace

OutputManager::OutputManager(const Particles& particles, const Input& input)
    : m_particles(particles), m_phi(input.modulationPotential()), m_id(input.id()), m_simname(input.simname()) {
  m_R = Utilities::LogAxis(input.ROutputMin(), input.ROutputMax(), input.ROutputSize());
}

double OutputManager::getFluxChargeGroup(int Z, double R) const {
  double value = 0.;
  for (const auto& particle : m_particles)
    if (particle.isChargeZ(Z)) value += particle.I_R_TOA(R, m_phi);
  return value;
}

double OutputManager::getFluxChargeIsotope(int Z, int A, double R) const {
  double value = 0.;
  for (const auto& particle : m_particles)
    if (particle.isChargeZ(Z) && particle.getPid().getA() == A) value += particle.I_R_TOA(R, m_phi);
  return value;
}

double OutputManager::getFluxChargeGroupEkn(int Z, double T) const {
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

void OutputManager::dumpSpectraRigidity() const {
  const std::string filename = "output/" + m_simname + "_spectra_R_" + std::to_string(m_id) + ".txt";
  std::ofstream out(filename);
  if (!out.is_open()) throw std::runtime_error("cannot open for writing: " + filename);
  LOGW << "writing rigidity spectra to " << filename;

  const double units = 1. / (CGS::GeV * pow2(CGS::meter) * CGS::sec);
  out << std::scientific;
  writeChargeGroupColumns(out, "R [GV]");
  for (const auto& R : m_R) {
    out << R / CGS::GeV << "\t";
    for (int Z = 1; Z <= 28; ++Z) out << getFluxChargeGroup(Z, R) / units << "\t";
    out << "\n";
  }
}

void OutputManager::dumpIsotopes() const {
  const std::string filename = "output/" + m_simname + "_isotopes_R_" + std::to_string(m_id) + ".txt";
  std::ofstream out(filename);
  if (!out.is_open()) throw std::runtime_error("cannot open for writing: " + filename);
  LOGW << "writing isotope spectra to " << filename;

  const double units = 1. / (CGS::GeV * pow2(CGS::meter) * CGS::sec);
  out << std::scientific;
  writeIsotopeColumns(out);
  for (const auto& R : m_R) {
    out << R / CGS::GeV << "\t";
    out << getFluxChargeIsotope(4, 7, R) / units << "\t";
    out << getFluxChargeIsotope(4, 9, R) / units << "\t";
    out << getFluxChargeIsotope(4, 10, R) / units << "\t";
    out << "\n";
  }
}

void OutputManager::dumpSpectraEkn() const {
  const std::string filename = "output/" + m_simname + "_spectra_Ekn_" + std::to_string(m_id) + ".txt";
  std::ofstream out(filename);
  if (!out.is_open()) throw std::runtime_error("cannot open for writing: " + filename);
  LOGW << "writing spectra to " << filename;

  const auto T = Utilities::LogAxis(0.3 * CGS::GeV, 1. * CGS::TeV, 4 * 32);
  const double units = 1. / (CGS::GeV * pow2(CGS::meter) * CGS::sec);
  out << std::scientific;
  writeChargeGroupColumns(out, "T [GeV/n]");
  for (const auto& T_i : T) {
    out << T_i / CGS::GeV << "\t";
    for (int Z = 1; Z <= 28; ++Z) out << getFluxChargeGroupEkn(Z, T_i) / units << "\t";
    out << "\n";
  }
}

}  // namespace CRAMS
