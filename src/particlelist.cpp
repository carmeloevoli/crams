#include "crams/particlelist.h"

#include <plog/Log.h>

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "crams/core/cgs.h"
#include "crams/utils/csvreader.h"
#include "crams/utils/utilities.h"

namespace {

const char kNucleilistFilename[] = "data/crams_nucleilist.csv";

struct AbundanceSetting {
  const char* key;
  int charge;
  double value;
  const char* label;
};

const AbundanceSetting kDefaultAbundances[] = {
    {"qh", 1, 5.06605e-02, "H"}, {"qhe", 2, 2.54369e-02, "He"},
    {"qli", 3, 0., "Li"},        {"qbe", 4, 0., "Be"},
    {"qb", 5, 0., "B"},          {"qc", 6, 3.98879e-03, "C"},
    {"qn", 7, 3.36117e-04, "N"}, {"qo", 8, 7.15129e-03, "O"},
    {"qf", 9, 0., "F"},          {"qne", 10, 1.34031e-03, "Ne"},
    {"qna", 11, 0.5e-4, "Na"},   {"qmg", 12, 2.38948e-03, "Mg"},
    {"qal", 13, 2.7e-4, "Al"},   {"qsi", 14, 2.77911e-03, "Si"},
    {"qp", 15, 1e-4, "P"},       {"qs", 16, 4.87000e-04, "S"},
    {"qcl", 17, 0., "Cl"},       {"qar", 18, 3e-4, "Ar"},
    {"qk", 19, 0., "K"},         {"qca", 20, 4e-4, "Ca"},
    {"qsc", 21, 0., "Sc"},       {"qti", 22, 0., "Ti"},
    {"qv", 23, 0., "V"},         {"qcr", 24, 2.5e-4, "Cr"},
    {"qmn", 25, 0., "Mn"},       {"qfe", 26, 6.80000e-03, "Fe"},
    {"qco", 27, 0., "Co"},       {"qni", 28, 4e-4, "Ni"},
};

struct SlopeSetting {
  const char* key;
  int charge;
  double value;
  const char* label;
};

const SlopeSetting kDefaultChargeSlopes[] = {
    {"hslope", 1, 4.37486, "H"},
    {"heslope", 2, 4.30995, "He"},
};

constexpr double kDefaultNucleiSlope = 4.32798;
constexpr int kDefaultNucleiMinCharge = 3;
constexpr size_t kNucleilistColumns = 5;

template <typename T>
T parseCsvValue(const std::vector<std::string>& row, size_t column, size_t rowIndex);

template <>
int parseCsvValue<int>(const std::vector<std::string>& row, size_t column, size_t rowIndex) {
  try {
    return std::stoi(row.at(column));
  } catch (const std::exception& e) {
    throw std::runtime_error("ParticleList: invalid integer in data row " + std::to_string(rowIndex) + ", column " +
                             std::to_string(column + 1) + ": " + e.what());
  }
}

template <>
double parseCsvValue<double>(const std::vector<std::string>& row, size_t column, size_t rowIndex) {
  try {
    return std::stod(row.at(column));
  } catch (const std::exception& e) {
    throw std::runtime_error("ParticleList: invalid floating-point value in data row " + std::to_string(rowIndex) +
                             ", column " + std::to_string(column + 1) + ": " + e.what());
  }
}

}  // namespace

namespace CRAMS {

ParticleList::ParticleList() {
  if (!Utilities::fileExists(kNucleilistFilename))
    throw std::runtime_error("ParticleList: nucleilist file not found: " + std::string(kNucleilistFilename));

  loadNucleilist(kNucleilistFilename);
  applyDefaultInjectionParameters();
}

ParticleList::ParticleList(const ParticleList& other) : m_list(other.m_list) { rebuildChargeIndex(); }

ParticleList& ParticleList::operator=(const ParticleList& other) {
  if (this != &other) {
    m_list = other.m_list;
    rebuildChargeIndex();
  }
  return *this;
}

ParticleList::ParticleList(ParticleList&& other) : m_list(std::move(other.m_list)) {
  rebuildChargeIndex();
  other.rebuildChargeIndex();
}

ParticleList& ParticleList::operator=(ParticleList&& other) {
  if (this != &other) {
    m_list = std::move(other.m_list);
    rebuildChargeIndex();
    other.rebuildChargeIndex();
  }
  return *this;
}

bool ParticleList::insert(const PID& key, const NucleusParameters& params) {
  const auto res = m_list.emplace(key, params);
  if (!res.second) {
    LOGD << "particle " << key << " already exists "
         << " with injection " << (res.first)->second;
    return false;
  }

  m_particlesByCharge[key.getZ()].push_back(res.first);
  return true;
}

void ParticleList::setAbundance(const PID& key, double value) {
  auto it = m_list.find(key);
  if (it != m_list.end()) {
    it->second.abundance = value;
    LOGD << "PID : " << key << " abundance modified to " << value;
  } else {
    LOGD << "PID : " << key << " not found in particle list";
  }
}

void ParticleList::setSlope(const PID& key, double value) {
  auto it = m_list.find(key);
  if (it != m_list.end()) {
    it->second.slope = value;
    LOGD << "PID : " << key << " slope modified to " << value;
  } else {
    LOGD << "PID : " << key << " not found in particle list";
  }
}

void ParticleList::setAbundanceChargeGroup(int charge, double abundance) {
  ensureChargeIndexFresh();
  const auto group = m_particlesByCharge.find(charge);
  if (group == m_particlesByCharge.end()) return;

  for (const auto& particle : group->second) {
    particle->second.abundance = abundance * particle->second.isotopicFractionISM;
  }
}

void ParticleList::setSlopeChargeGroup(int charge, double slope) {
  ensureChargeIndexFresh();
  const auto group = m_particlesByCharge.find(charge);
  if (group == m_particlesByCharge.end()) return;

  for (const auto& particle : group->second) {
    particle->second.slope = slope;
  }
}

void ParticleList::setSlopeNuclei(int minCharge, double slope) {
  ensureChargeIndexFresh();
  for (auto group = m_particlesByCharge.lower_bound(minCharge); group != m_particlesByCharge.end(); ++group) {
    for (const auto& particle : group->second) {
      particle->second.slope = slope;
    }
  }
}

void ParticleList::setParam(const std::string& key, double value) {
  const auto simpleKey = Utilities::simplifyKey(key);

  for (const auto& setting : kDefaultAbundances) {
    if (simpleKey == setting.key) {
      setAbundanceChargeGroup(setting.charge, value);
      LOGD << "changed " << setting.label << " abundance to " << value;
      return;
    }
  }

  for (const auto& setting : kDefaultChargeSlopes) {
    if (simpleKey == setting.key) {
      setSlopeChargeGroup(setting.charge, value);
      LOGD << "changed " << setting.label << " slope to " << value;
      return;
    }
  }

  if (simpleKey == "slope") {
    setSlopeNuclei(kDefaultNucleiMinCharge, value);
    LOGD << "changed nuclei slope to " << value;
    return;
  }

  LOGD << "ignored unknown particle parameter " << key;
}

void ParticleList::readParamsFromFile(const std::string& filename) {
  std::ifstream infile(filename);
  if (!infile.is_open()) throw std::runtime_error("ParticleList: cannot open parameter file '" + filename + "'");

  std::string line;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);
    std::string key;
    double value;
    if (!(iss >> key >> value)) continue;
    setParam(key, value);
  }
}

void ParticleList::loadNucleilist(const std::string& filename) {
  CSVReader reader(filename);
  const auto nucleilist = reader.getData();

  m_list.clear();
  m_particlesByCharge.clear();
  m_rebuildChargeIndexBeforeUpdate = false;

  for (size_t i = 0; i < nucleilist.size(); ++i) {
    const auto& row = nucleilist[i];
    const size_t rowIndex = i + 1;
    if (row.size() < kNucleilistColumns) {
      throw std::runtime_error("ParticleList: expected at least " + std::to_string(kNucleilistColumns) +
                               " columns in data row " + std::to_string(rowIndex));
    }

    const int Z = parseCsvValue<int>(row, 0, rowIndex);
    const int A = parseCsvValue<int>(row, 1, rowIndex);
    const bool isTertiary = parseCsvValue<int>(row, 2, rowIndex) != 0;
    const double decayHalfLife = parseCsvValue<double>(row, 3, rowIndex) * CGS::Myr;
    const double isotopicFractionISM = 0.01 * parseCsvValue<double>(row, 4, rowIndex);

    // Injection abundance and slope are set by applyDefaultInjectionParameters()
    // and the .ini, so they start at zero here.
    const auto pid = PID{Z, A, isTertiary};
    const auto params = NucleusParameters{0., 0., isotopicFractionISM, decayHalfLife, decayHalfLife < 0., false};
    insert(pid, params);
  }
}

void ParticleList::applyDefaultInjectionParameters() {
  for (const auto& setting : kDefaultAbundances) {
    setAbundanceChargeGroup(setting.charge, setting.value);
  }

  for (const auto& setting : kDefaultChargeSlopes) {
    setSlopeChargeGroup(setting.charge, setting.value);
  }

  setSlopeNuclei(kDefaultNucleiMinCharge, kDefaultNucleiSlope);
}

void ParticleList::ensureChargeIndexFresh() {
  if (m_rebuildChargeIndexBeforeUpdate) rebuildChargeIndex();
}

void ParticleList::rebuildChargeIndex() {
  m_particlesByCharge.clear();
  for (auto particle = m_list.begin(); particle != m_list.end(); ++particle) {
    m_particlesByCharge[particle->first.getZ()].push_back(particle);
  }
}

void ParticleList::print() const {
  LOGI << "Particle list contains " << m_list.size() << " nuclei.";
  for (const auto& particle : m_list) LOGD << "found nucleus " << particle.first << " with params " << particle.second;
}

}  // namespace CRAMS
