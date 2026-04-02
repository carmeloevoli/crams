#include "chi2.h"

#include <plog/Log.h>

#define max_num_of_char_in_a_line 512
#define num_of_header_lines 7

namespace CRAMS {

Chi2::Chi2() {}

Chi2::Chi2(const std::string& name, const Particles& particles, const double& phi)
    : m_name(name), m_particles(particles), m_phi(phi) {}

void Chi2::readKISSfile(const std::string& filename, const double& xunits, const double& yunits) {
  LOGI << "reading data from " << filename << "... ";
  std::ifstream file_to_read(filename.c_str());
  if (file_to_read.is_open()) {
    for (int i = 0; i < num_of_header_lines; ++i) file_to_read.ignore(max_num_of_char_in_a_line, '\n');
    file_to_read.ignore(max_num_of_char_in_a_line, '\n');
    double values[6];
    while (!file_to_read.eof()) {
      file_to_read >> values[0] >> values[1] >> values[2] >> values[3] >> values[4] >> values[5];
      dataPoint point;
      point.R = values[0] * xunits;
      point.I = values[1] * yunits;
      point.I_err = std::make_pair<double, double>(values[2] * yunits, values[3] * yunits);
      if (file_to_read.good()) m_data.push_back(point);
    }
  } else {
    throw std::runtime_error("data file cannot be open!");
  }
  LOGI << " with size : " << m_data.size();
  file_to_read.close();
}

void Chi2::readCRDBfile(const std::string& filename, const double& xunits, const double& yunits) {
  LOGI << "reading data from " << filename << "... ";
  std::ifstream file_to_read(filename.c_str());
  if (file_to_read.is_open()) {
    for (int i = 0; i < num_of_header_lines; ++i) file_to_read.ignore(max_num_of_char_in_a_line, '\n');
    file_to_read.ignore(max_num_of_char_in_a_line, '\n');
    double values[6];
    while (!file_to_read.eof()) {
      file_to_read >> values[0] >> values[1] >> values[2] >> values[3] >> values[4] >> values[5];
      dataPoint point;
      point.R = values[0] * xunits;
      point.I = values[1] * yunits;
      point.I_err = std::make_pair<double, double>(std::sqrt(values[2]*values[2]+values[4]*values[4]) * yunits, std::sqrt(values[3]*values[3]+values[5]*values[5]) * yunits);
      if (file_to_read.good()) m_data.push_back(point);
    }
  } else {
    throw std::runtime_error("data file cannot be open!");
  }
  LOGI << " with size : " << m_data.size();
  file_to_read.close();
}

double Chi2::computeChi2(const double& R_min, const double& R_max) const {
  // If using kinetic energy, delegate to the kinetic energy computation
  if (m_energyType == EnergyType::KINETIC_ENERGY_PER_NUCLEON || 
      m_energyType == EnergyType::KINETIC_ENERGY_TOTAL) {
    return computeChi2KineticEnergy(R_min, R_max);
  }
  
  // Otherwise use rigidity-based computation
  double chi2 = 0.0;
  size_t ndata = 0;
  for (auto idata = m_data.begin(); idata != m_data.end(); idata++) {
    if (idata->R > R_min && idata->R < R_max) {
      const double I_R_TOA = getModel(idata->R, m_phi);
      double delta_chi2 = pow2(I_R_TOA - idata->I);
      delta_chi2 /= (I_R_TOA < idata->I) ? pow2(idata->I_err.first) : pow2(idata->I_err.second);
      chi2 += delta_chi2;
      ndata++;
    }
  }
  return chi2 / (double)ndata;
}

double Chi2::computeChi2KineticEnergy(const double& E_min, const double& E_max) const {
  double chi2 = 0.0;
  size_t ndata = 0;
  if (m_energyType == EnergyType::KINETIC_ENERGY_PER_NUCLEON) {
    for (auto idata = m_data.begin(); idata != m_data.end(); idata++) {
      if (idata->R > E_min && idata->R < E_max) {
        // For kinetic energy data, R field stores the energy value (Ekn or Ek)
        const double model_intensity = getModel(idata->R, m_phi);
        double delta_chi2 = pow2(model_intensity - idata->I);
        delta_chi2 /= (model_intensity < idata->I) ? pow2(idata->I_err.first) : pow2(idata->I_err.second);
        chi2 += delta_chi2;
        ndata++;
      }
    }
  } else if (m_energyType == EnergyType::KINETIC_ENERGY_TOTAL) {
    for (auto idata = m_data.begin(); idata != m_data.end(); idata++) {
        // For total kinetic energy data, R field stores the total energy value (Ek)
        const double model_intensity = getModel(idata->R, m_phi);
        double delta_chi2 = pow2(model_intensity - idata->I);
        delta_chi2 /= (model_intensity < idata->I) ? pow2(idata->I_err.first) : pow2(idata->I_err.second);
        chi2 += delta_chi2;
        ndata++;
    }
  } else {
    LOGW << "Energy type is not kinetic energy, cannot compute chi2 for kinetic energy data!";
    return 0.0;
  }
  
  return (ndata > 0) ? chi2 / (double)ndata : 0.0;
}

double Chi2IH::getModel(const double& R, const double& phi) const {
  double value = (itH1.first) ? itH1.second->I_R_TOA(R, phi) : 0.;
  value += (itH2.first) ? itH2.second->I_R_TOA(R, phi) : 0.;
  value += (itH1_ter.first) ? itH1_ter.second->I_R_TOA(R, phi) : 0.;
  return value;
}

double Chi2IHe::getModel(const double& R, const double& phi) const {
  double value = (itHe3.first) ? itHe3.second->I_R_TOA(R, phi) : 0.;
  value += (itHe4.first) ? itHe4.second->I_R_TOA(R, phi) : 0.;
  return value;
}

double Chi2BC::getModel(const double& R, const double& phi) const {
  double B = (itB10.first) ? itB10.second->I_R_TOA(R, phi) : 0.;
  B += (itB11.first) ? itB11.second->I_R_TOA(R, phi) : 0.;
  double C = (itC12.first) ? itC12.second->I_R_TOA(R, phi) : 0.;
  C += (itC13.first) ? itC13.second->I_R_TOA(R, phi) : 0.;
  C += (itC14.first) ? itC14.second->I_R_TOA(R, phi) : 0.;
  return B / C;
}

// ============================================================================
// Chi2Composite implementation
// ============================================================================

double Chi2Composite::evaluateGroup(const std::vector<IsotopeGroup>& groups,
                                     const double& energyVal,
                                     const double& phi) const {
  double result = 0.0;
  
  for (const auto& group : groups) {
    double groupValue = 0.0;
    
    for (const auto& pid : group.pids) {
      if (m_cachedIterators.find(pid) == m_cachedIterators.end()) {
        continue;  // PID not found, skip
      }
      
      const auto& it = m_cachedIterators.at(pid);
      if (!it.first) continue;  // Particle not present
      
      double intensity = 0.0;
      
      if (m_energyType == EnergyType::RIGIDITY) {
        // Use rigidity for AMS-02 style data
        intensity = it.second->I_R_TOA(energyVal, phi);
      } else if (m_energyType == EnergyType::KINETIC_ENERGY_PER_NUCLEON) {
        // Use kinetic energy per nucleon for DAMPE flux data
        // energyVal is already Ekn in GeV, convert to rigidity equivalent
        double Ekn = energyVal;  // Kinetic energy per nucleon
        intensity = it.second->I_T_TOA(Ekn, phi);
      } else if (m_energyType == EnergyType::KINETIC_ENERGY_TOTAL) {
        // Use total kinetic energy for DAMPE ratio data
        double Ek = energyVal;  // Total kinetic energy
        int A = pid.getA();  // Mass number
        if (A > 0) {
          double Ekn = Ek / A;  // Convert total energy to per-nucleon
          intensity = it.second->I_T_TOA(Ekn, phi)/A;
        }
      }
      
      groupValue += intensity;
    }
    
    // Apply group scaling
    groupValue *= group.scaling;
    result += groupValue;
  }
  
  return result;
}

double Chi2Composite::getModel(const double& energyVal, const double& phi) const {
  double numerator = evaluateGroup(m_numeratorGroups, energyVal, phi);
  
  if (m_mode == CompositionMode::SUM) {
    return numerator;
  } else if (m_mode == CompositionMode::RATIO) {
    double denominator = evaluateGroup(m_denominatorGroups, energyVal, phi);
    if (denominator == 0.0) return 0.0;
    return numerator / denominator;
  } else if (m_mode == CompositionMode::SCALED) {
    // For SCALED mode, apply energy-dependent scaling
    double scaling = std::pow(energyVal, 2.7);  // Example: can be customized
    return numerator * scaling;
  }
  
  return numerator;
}

// ============================================================================
// Factory function implementations
// ============================================================================

std::string getElementSymbol(int Z) {
  switch(Z) {
    case 1: return "H";
    case 2: return "He";
    case 3: return "Li";
    case 4: return "Be";
    case 5: return "B";
    case 6: return "C";
    case 7: return "N";
    case 8: return "O";
    case 9: return "F";
    case 10: return "Ne";
    case 11: return "Na";
    case 12: return "Mg";
    case 13: return "Al";
    case 14: return "Si";
    case 15: return "P";
    case 16: return "S";
    case 17: return "Cl";
    case 18: return "Ar";
    case 19: return "K";
    case 20: return "Ca";
    case 21: return "Sc";
    case 22: return "Ti";
    case 23: return "V";
    case 24: return "Cr";
    case 25: return "Mn";
    case 26: return "Fe";
    case 27: return "Co";
    case 28: return "Ni";
    default: return "";
  }
}

// Manual version of makeElement (isotope vector)
std::unique_ptr<Chi2> get_Chi2(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    const std::vector<PID>& isotopes,
    const std::string& dataFile,
    const double& xunits,
    const double& yunits) {
  std::vector<IsotopeGroup> numerator;
  numerator.push_back({isotopes, 1.0, false});
  
  return std::make_unique<Chi2Composite>(
      name, particles, phi, 
      CompositionMode::SUM,
      EnergyType::RIGIDITY,
      dataFile, xunits, yunits,
      numerator);
}

// Automatic version of makeElement (Z atomic number)
std::unique_ptr<Chi2> get_Chi2(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    int Z,
    const std::string& dataFile,
    const double& xunits,
    const double& yunits) {
  std::vector<PID> isotopes = findIsotopesOfElement(Z);
  
  if (isotopes.empty()) {
    throw std::runtime_error("No isotopes found for element with Z=" + std::to_string(Z));
  }
  
  // Auto-find data file if not provided
  std::string actualDataFile = dataFile;
  if (actualDataFile.empty()) {
    std::string elementSymbol = getElementSymbol(Z);
    if (!elementSymbol.empty()) {
      actualDataFile = findDataFile(elementSymbol);
    }
  }
  
  LOGI << "Creating element model '" << name << "' for Z=" << Z << " with " << isotopes.size() << " isotopes";
  
  return get_Chi2(name, particles, phi, isotopes, actualDataFile, xunits, yunits);
}

// ============================================================================
// Automatic isotope finding implementations
// ============================================================================

std::vector<PID> findIsotopesOfElement(int Z) {
  std::vector<PID> isotopes;
  
  // This is the definitive list of all supported isotopes from pid.h
  // organized by atomic number (Z)
  switch(Z) {
    case 1:  // Hydrogen
      isotopes = {H1, H2, H1_ter};
      break;
    case 2:  // Helium
      isotopes = {He3, He4};
      break;
    case 3:  // Lithium
      isotopes = {Li6, Li7};
      break;
    case 4:  // Beryllium
      isotopes = {Be7, Be9, Be10};
      break;
    case 5:  // Boron
      isotopes = {B10, B11};
      break;
    case 6:  // Carbon
      isotopes = {C12, C13, C14};
      break;
    case 7:  // Nitrogen
      isotopes = {N14, N15};
      break;
    case 8:  // Oxygen
      isotopes = {O16, O17, O18};
      break;
    case 9:  // Fluorine
      isotopes = {F19};
      break;
    case 10:  // Neon
      isotopes = {Ne20, Ne21, Ne22};
      break;
    case 11:  // Sodium
      isotopes = {Na22, Na23};
      break;
    case 12:  // Magnesium
      isotopes = {Mg24, Mg25, Mg26};
      break;
    case 13:  // Aluminum
      isotopes = {Al26, Al27};
      break;
    case 14:  // Silicon
      isotopes = {Si28, Si29, Si30, Si32};
      break;
    case 15:  // Phosphorus
      isotopes = {P31, P32, P33};
      break;
    case 16:  // Sulfur
      isotopes = {S32, S33, S34, S36};
      break;
    case 17:  // Chlorine
      isotopes = {Cl35, Cl36, Cl37};
      break;
    case 18:  // Argon
      isotopes = {Ar36, Ar37, Ar38, Ar40};
      break;
    case 19:  // Potassium
      isotopes = {K39, K40, K41};
      break;
    case 20:  // Calcium
      isotopes = {Ca40, Ca41, Ca42, Ca43, Ca44, Ca46, Ca48};
      break;
    case 21:  // Scandium
      isotopes = {Sc45};
      break;
    case 22:  // Titanium
      isotopes = {Ti44, Ti46, Ti47, Ti48, Ti49, Ti50};
      break;
    case 23:  // Vanadium
      isotopes = {V49, V50, V51};
      break;
    case 24:  // Chromium
      isotopes = {Cr48, Cr50, Cr51, Cr52, Cr53, Cr54};
      break;
    case 25:  // Manganese
      isotopes = {Mn53, Mn54, Mn55};
      break;
    case 26:  // Iron
      isotopes = {Fe54, Fe55, Fe56, Fe57, Fe58, Fe60};
      break;
    case 27:  // Cobalt
      isotopes = {Co57, Co59};
      break;
    case 28:  // Nickel
      isotopes = {Ni56, Ni58, Ni59, Ni60, Ni61, Ni62, Ni64};
      break;
    default:
      LOGW << "Element with Z=" << Z << " not found in predefined isotope list";
      break;
  }
  
  return isotopes;
}

// Manual version of get_Chi2_ratio (isotope vectors)
std::unique_ptr<Chi2> get_Chi2_ratio(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    const std::vector<PID>& numeratorIsotopes,
    const std::vector<PID>& denominatorIsotopes,
    const std::string& dataFile,
    const double& xunits,
    const double& yunits) {
  std::vector<IsotopeGroup> numerator;
  numerator.push_back({numeratorIsotopes, 1.0, false});
  
  std::vector<IsotopeGroup> denominator;
  denominator.push_back({denominatorIsotopes, 1.0, false});
  
  return std::make_unique<Chi2Composite>(
      name, particles, phi,
      CompositionMode::RATIO,
      EnergyType::RIGIDITY,
      dataFile, xunits, yunits,
      numerator, denominator);
}

// Automatic version of get_Chi2_ratio (Z atomic numbers)
std::unique_ptr<Chi2> get_Chi2_ratio(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    int numeratorZ,
    int denominatorZ,
    const std::string& dataFile,
    const double& xunits,
    const double& yunits) {
  std::vector<PID> numeratorIsotopes = findIsotopesOfElement(numeratorZ);
  std::vector<PID> denominatorIsotopes = findIsotopesOfElement(denominatorZ);
  
  if (numeratorIsotopes.empty()) {
    throw std::runtime_error("No isotopes found for numerator element with Z=" + std::to_string(numeratorZ));
  }
  if (denominatorIsotopes.empty()) {
    throw std::runtime_error("No isotopes found for denominator element with Z=" + std::to_string(denominatorZ));
  }
  
  // Auto-find data file if not provided
  std::string actualDataFile = dataFile;
  if (actualDataFile.empty()) {
    std::string numeratorSymbol = getElementSymbol(numeratorZ);
    std::string denominatorSymbol = getElementSymbol(denominatorZ);
    if (!numeratorSymbol.empty() && !denominatorSymbol.empty()) {
      actualDataFile = findDataFile(numeratorSymbol + "_" + denominatorSymbol);
    }
  }
  
  LOGI << "Creating ratio model '" << name << "' with Z=" << numeratorZ << " (" << numeratorIsotopes.size() 
       << " isotopes) / Z=" << denominatorZ << " (" << denominatorIsotopes.size() << " isotopes)";
  
  return get_Chi2_ratio(name, particles, phi, numeratorIsotopes, denominatorIsotopes, actualDataFile, xunits, yunits);
}

// Create Chi2 for DAMPE total flux data (as function of total kinetic energy Ek in GeV)
std::unique_ptr<Chi2> get_Chi2_DAMPE_flux(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    int Z,
    const std::string& dataFile,
    const double& xunits,
    const double& yunits) {
  std::vector<PID> isotopes = findIsotopesOfElement(Z);
  
  if (isotopes.empty()) {
    throw std::runtime_error("No isotopes found for element with Z=" + std::to_string(Z));
  }
  
  // Auto-find data file if not provided
  std::string actualDataFile = dataFile;
  if (actualDataFile.empty()) {
    std::string elementSymbol = getElementSymbol(Z);
    if (!elementSymbol.empty()) {
      // Look for DAMPE_Element_Ekn.txt or similar
      std::vector<std::string> patterns = {
          "data/DAMPE_" + elementSymbol + "_Ekn.txt",
          "data/DAMPE_" + elementSymbol + "_Ek.txt",
          "DAMPE_" + elementSymbol + "_Ekn.txt",
          "DAMPE_" + elementSymbol + "_Ek.txt"
      };
      
      for (const auto& pattern : patterns) {
        if (fileExists(pattern)) {
          actualDataFile = pattern;
          break;
        }
      }
    }
  }
  
  LOGI << "Creating DAMPE flux model '" << name << "' for Z=" << Z << " with " << isotopes.size() << " isotopes";
  
  std::vector<IsotopeGroup> numerator;
  numerator.push_back({isotopes, 1.0, false});
  
  return std::make_unique<Chi2Composite>(
      name, particles, phi,
      CompositionMode::SUM,
      EnergyType::KINETIC_ENERGY_TOTAL,
      actualDataFile, xunits, yunits,
      numerator);
}

// Create Chi2 for DAMPE ratio data (as function of kinetic energy per nucleon Ekn in GeV)
std::unique_ptr<Chi2> get_Chi2_DAMPE_ratio(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    int numeratorZ,
    int denominatorZ,
    const std::string& dataFile,
    const double& xunits,
    const double& yunits) {
  std::vector<PID> numeratorIsotopes = findIsotopesOfElement(numeratorZ);
  std::vector<PID> denominatorIsotopes = findIsotopesOfElement(denominatorZ);
  
  if (numeratorIsotopes.empty()) {
    throw std::runtime_error("No isotopes found for numerator element with Z=" + std::to_string(numeratorZ));
  }
  if (denominatorIsotopes.empty()) {
    throw std::runtime_error("No isotopes found for denominator element with Z=" + std::to_string(denominatorZ));
  }
  
  // Auto-find data file if not provided
  std::string actualDataFile = dataFile;
  if (actualDataFile.empty()) {
    std::string numeratorSymbol = getElementSymbol(numeratorZ);
    std::string denominatorSymbol = getElementSymbol(denominatorZ);
    if (!numeratorSymbol.empty() && !denominatorSymbol.empty()) {
      // Look for DAMPE_Numerator_Denominator_Ekn.txt or similar
      std::vector<std::string> patterns = {
          "data/DAMPE_" + numeratorSymbol + "_" + denominatorSymbol + "_Ekn.txt",
          "DAMPE_" + numeratorSymbol + "_" + denominatorSymbol + "_Ekn.txt"
      };
      
      for (const auto& pattern : patterns) {
        if (fileExists(pattern)) {
          actualDataFile = pattern;
          break;
        }
      }
    }
  }
  
  LOGI << "Creating DAMPE ratio model '" << name << "' with Z=" << numeratorZ << " (" << numeratorIsotopes.size() 
       << " isotopes) / Z=" << denominatorZ << " (" << denominatorIsotopes.size() << " isotopes)";
  
  std::vector<IsotopeGroup> numerator;
  numerator.push_back({numeratorIsotopes, 1.0, false});
  
  std::vector<IsotopeGroup> denominator;
  denominator.push_back({denominatorIsotopes, 1.0, false});
  
  return std::make_unique<Chi2Composite>(
      name, particles, phi,
      CompositionMode::RATIO,
      EnergyType::KINETIC_ENERGY_PER_NUCLEON,
      actualDataFile, xunits, yunits,
      numerator, denominator);
}

}  // namespace CRAMS

// double Chi2_B::get_model(const double& R, const double& phi) const {
// 	double value = (ptr_B10.isPresent) ? ptr_B10.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_B11.isPresent) ? ptr_B11.it->I_R_TOA(R, phi) : 0.;
// 	return value;
// }

// double Chi2_C::get_model(const double& R, const double& phi) const {
// 	double value = (ptr_C12.isPresent) ? ptr_C12.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_C13.isPresent) ? ptr_C13.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_C14.isPresent) ? ptr_C14.it->I_R_TOA(R, phi) : 0.;
// 	return value;
// }

// double Chi2_N::get_model(const double& R, const double& phi) const {
// 	double value = (ptr_N14.isPresent) ? ptr_N14.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_N15.isPresent) ? ptr_N15.it->I_R_TOA(R, phi) : 0.;
// 	return value;
// }

// double Chi2_O::get_model(const double& R, const double& phi) const {
// 	double value = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
// 	value += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
// 	return value;
// }

// double Chi2_HeO::get_model(const double& R, const double& phi) const {
// 	double O = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
// 	double He = (ptr_He3.isPresent) ? ptr_He3.it->I_R_TOA(R, phi) : 0.;
// 	He += (ptr_He4.isPresent) ? ptr_He4.it->I_R_TOA(R, phi) : 0.;
// 	return He / O;
// }

// double Chi2_BeB::get_model(const double& R, const double& phi) const {
// 	double B = (ptr_B10.isPresent) ? ptr_B10.it->I_R_TOA(R, phi) : 0.;
// 	B += (ptr_B11.isPresent) ? ptr_B11.it->I_R_TOA(R, phi) : 0.;
// 	double Be = (ptr_Be7.isPresent) ? ptr_Be7.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be9.isPresent) ? ptr_Be9.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be10.isPresent) ? ptr_Be10.it->I_R_TOA(R, phi) : 0.;
// 	return Be / B;
// }

// double Chi2_BeB_statsonly::get_model(const double& R, const double& phi) const {
// 	double B = (ptr_B10.isPresent) ? ptr_B10.it->I_R_TOA(R, phi) : 0.;
// 	B += (ptr_B11.isPresent) ? ptr_B11.it->I_R_TOA(R, phi) : 0.;
// 	double Be = (ptr_Be7.isPresent) ? ptr_Be7.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be9.isPresent) ? ptr_Be9.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be10.isPresent) ? ptr_Be10.it->I_R_TOA(R, phi) : 0.;
// 	return Be / B;
// }

// double Chi2_BeC::get_model(const double& R, const double& phi) const {
// 	double C = (ptr_C12.isPresent) ? ptr_C12.it->I_R_TOA(R, phi) : 0.;
// 	C += (ptr_C13.isPresent) ? ptr_C13.it->I_R_TOA(R, phi) : 0.;
// 	C += (ptr_C14.isPresent) ? ptr_C14.it->I_R_TOA(R, phi) : 0.;
// 	double Be = (ptr_Be7.isPresent) ? ptr_Be7.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be9.isPresent) ? ptr_Be9.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be10.isPresent) ? ptr_Be10.it->I_R_TOA(R, phi) : 0.;
// 	return Be / C;
// }

// double Chi2_BeO::get_model(const double& R, const double& phi) const {
// 	double O = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
// 	double Be = (ptr_Be7.isPresent) ? ptr_Be7.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be9.isPresent) ? ptr_Be9.it->I_R_TOA(R, phi) : 0.;
// 	Be += (ptr_Be10.isPresent) ? ptr_Be10.it->I_R_TOA(R, phi) : 0.;
// 	return Be / O;
// }

// double Chi2_BO::get_model(const double& R, const double& phi) const {
// 	double B = (ptr_B10.isPresent) ? ptr_B10.it->I_R_TOA(R, phi) : 0.;
// 	B += (ptr_B11.isPresent) ? ptr_B11.it->I_R_TOA(R, phi) : 0.;
// 	double O = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
// 	return B / O;
// }

// double Chi2_CO::get_model(const double& R, const double& phi) const {
// 	double O = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
// 	O += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
// 	double C = (ptr_C12.isPresent) ? ptr_C12.it->I_R_TOA(R, phi) : 0.;
// 	C += (ptr_C13.isPresent) ? ptr_C13.it->I_R_TOA(R, phi) : 0.;
// 	C += (ptr_C14.isPresent) ? ptr_C14.it->I_R_TOA(R, phi) : 0.;
// 	return C / O;
// }
