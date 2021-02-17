#include "particlelist.h"

#include <plog/Log.h>

#include <sstream>

#include "cgs.h"
#include "csvreader.h"
#include "utilities.h"

namespace CRAMS {

ParticleList::ParticleList() {
  if (Utilities::fileExists(nucleilistFilename))
    loadNucleilist(nucleilistFilename);
  else
    throw std::runtime_error("nuclelist file not found");
}

ParticleList::~ParticleList() {
  m_list.clear();
  LOGD << "released memory from ParticleList";
}

bool ParticleList::insert(const PID& key, const NucleusParameters& params) {
  auto res = m_list.insert(std::make_pair(key, params));
#ifdef DEBUG
  if (!res.second) {
    std::cout << "particle " << key << " already exists "
              << " with injection " << (res.first)->second << std::endl;
  } else {
    // spdlog::debug("created particle {}", 1);  // << key << " with parameters " << params << std::endl;
  }
#endif
  return res.second;
}

void ParticleList::setAbundance(const PID& key, const double& value) {
  auto it = m_list.find(key);
  if (it != m_list.end()) {
    it->second.abundance = value;
#ifdef DEBUG
    std::cout << "PID : " << key << " abundance modified to " << value << ".\n";
#endif
  } else {
#ifdef DEBUG
    std::cout << "PID : " << key << " not found in particle list.\n";
#endif
  }
}

void ParticleList::setSlope(const PID& key, const double& value) {
  auto it = m_list.find(key);
  if (it != m_list.end()) {
    it->second.slope = value;
#ifdef DEBUG
    std::cout << "PID : " << key << " slope modified to " << value << ".\n";
#endif
  } else {
#ifdef DEBUG
    std::cout << "PID : " << key << " not found in particle list.\n";
#endif
  }
}

// void ParticleList::set_from_file(const std::string& filename) {
//   const std::ifstream infile(filename.c_str());
//   std::string line;
//   while (std::getline(infile, line)) {
//     std::istringstream iss(line);
//     std::string key;
//     double value;
//     if (!(iss >> key >> value)) {
//       break;
//     }  // error
//     if (key == "q_H")
//       set_abundance(H1, value);
//     else if (key == "q_He")
//       set_abundance(He4, value);
//     else if (key == "q_Li")
//       set_abundance(Li6, value);
//     else if (key == "q_Be")
//       set_abundance(Be9, value);
//     else if (key == "q_B")
//       set_abundance(B10, value);
//     else if (key == "q_C")
//       set_abundance(C12, value);
//     else if (key == "q_N")
//       set_abundance(N14, value);
//     else if (key == "q_O")
//       set_abundance(O16, value);
//   }
// }

void ParticleList::setAbundanceChargeGroup(const int& charge, const double& abundance) {
  for (auto& p : m_list) {
    if (p.first.getZ() == charge) p.second.abundance = abundance * p.second.isotopicFractionISM;
  }
}

void ParticleList::setSlopeChargeGroup(const int& charge, const double& slope) {
  for (auto& p : m_list) {
    if (p.first.getZ() == charge) p.second.slope = slope;
  }
}

void ParticleList::setSlopeNuclei(const int& minCharge, const double& slope) {
  for (auto& p : m_list) {
    if (p.first.getZ() >= minCharge) p.second.slope = slope;
  }
}

void ParticleList::readParamsFromFile(const std::string& filename) {
  setAbundanceChargeGroup(1, 5.04913e-02);   // q_H
  setAbundanceChargeGroup(2, 2.51187e-02);   // q_He
  setAbundanceChargeGroup(3, 0.);            // q_Li
  setAbundanceChargeGroup(4, 0.);            // q_Be
  setAbundanceChargeGroup(5, 0.);            // q_B
  setAbundanceChargeGroup(6, 3.98879e-03);   // q_C
  setAbundanceChargeGroup(7, 3.36117e-04);   // q_N
  setAbundanceChargeGroup(8, 7.15129e-03);   // q_O
  setAbundanceChargeGroup(9, 0.);            // q_F
  setAbundanceChargeGroup(10, 1.34031e-03);  // q_Ne
  setAbundanceChargeGroup(11, .5e-4);        // q_Na
  setAbundanceChargeGroup(12, 2.38948e-03);  // q_Mg
  setAbundanceChargeGroup(13, 2.7e-4);       // q_Al
  setAbundanceChargeGroup(14, 2.77911e-03);  // q_Si
  setAbundanceChargeGroup(15, 1e-4);         // q_P
  setAbundanceChargeGroup(16, 4.80499e-04);  // q_S
  setAbundanceChargeGroup(17, 0.);           // q_Cl
  setAbundanceChargeGroup(18, 3e-4);         // q_Ar
  setAbundanceChargeGroup(19, 0.);           // q_K
  setAbundanceChargeGroup(20, 4e-4);         // q_Ca
  setAbundanceChargeGroup(21, 0.);           // q_Sc
  setAbundanceChargeGroup(22, 0.);           // q_Ti
  setAbundanceChargeGroup(23, 0.);           // q_V
  setAbundanceChargeGroup(24, 2.5e-4);       // q_Cr
  setAbundanceChargeGroup(25, 0.);           // q_Mn
  setAbundanceChargeGroup(26, 6.9e-3);       // q_Fe
  setAbundanceChargeGroup(27, 0.);           // q_Co
  setAbundanceChargeGroup(28, 4e-4);         // q_Ni
  setSlopeChargeGroup(1, 4.37326e+00);
  setSlopeChargeGroup(2, 4.30578e+00);
  setSlopeNuclei(3, 4.33134);
}

void ParticleList::loadNucleilist(const std::string& filename) {
  CSVReader reader(filename);
  auto nucleilist = reader.getData();
  for (const auto& vec : nucleilist) {
    const int Z = atoi(vec[0].c_str());
    const int A = atoi(vec[1].c_str());
    const bool isTertiary = atoi(vec[2].c_str());
    const double decayTime = atof(vec[3].c_str()) * CGS::Myr;
    const double isotopicFractionISM = 0.01 * atof(vec[4].c_str());
    const double injectionAb = atof(vec[5].c_str());
    const double injectionSlope = atof(vec[6].c_str());
    auto pid = PID{Z, A, isTertiary};
    bool isStable = (decayTime < 0.) ? true : false;
    auto params = NucleusParameters{injectionAb, injectionSlope, isotopicFractionISM, decayTime, isStable, false};
    insert(pid, params);
  }
}

void ParticleList::print() const {
  LOGI << "Particle list contains " << m_list.size() << " nuclei.";
  for (auto& particle : m_list) LOGD << "found nucleus " << particle.first << " with params " << particle.second;
}

}  // namespace CRAMS