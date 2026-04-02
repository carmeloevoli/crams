#ifndef INCLUDE_CHI2_H_
#define INCLUDE_CHI2_H_

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "particle.h"
#include "pid.h"

namespace CRAMS {

struct dataPoint {
  double R;
  double I;
  std::pair<double, double> I_err;
  int A = 0;  // Mass number (for per-nucleon kinetic energy data)
};

// Operation modes for Chi2Composite
enum class CompositionMode {
  SUM,     // Add all isotope intensities
  RATIO,   // Divide numerator group by denominator group
  SCALED   // Scale result by energy-dependent factor
};

// Energy type: Rigidity (R) or Kinetic energy (T)
enum class EnergyType {
  RIGIDITY,              // I_R_TOA(R, phi) - rigidity in MeV
  KINETIC_ENERGY_PER_NUCLEON,  // I_T_TOA(T/A, phi) - kinetic energy per nucleon in GeV (DAMPE flux files)
  KINETIC_ENERGY_TOTAL   // I_T_TOA(T, phi) - total kinetic energy in GeV (DAMPE ratio files)
};

class Chi2 {
 public:
  Chi2();
  Chi2(const std::string& name, const Particles& particles, const double& phi);
  virtual ~Chi2() = default;
  void setPhi(const double& modulationPotential) { m_phi = modulationPotential; }
  inline double getPhi() const { return m_phi; }
  inline double getChi2() const { return m_chi2; }
  inline std::string getName() const { return m_name; }
  double computeChi2(const double& R_min, const double& R_max = 5. * CGS::TeV) const;
  double computeChi2KineticEnergy(const double& E_min, const double& E_max) const;

  itParticle findParticle(const PID& pid) {
    auto it = find(m_particles.begin(), m_particles.end(), Particle(pid));
    bool isPresent = !(it == m_particles.end());
    return itParticle(isPresent, it);
  }

 protected:
  virtual double getModel(const double& R, const double& phi) const { return 0; }
  void readKISSfile(const std::string& filename, const double& xunits, const double& yunits = 1.0);
  void readCRDBfile(const std::string& filename, const double& xunits, const double& yunits = 1.0);

 protected:
  const double xunits = CGS::GeV;
  std::string m_name;
  double m_chi2 = 0;
  double m_phi = 0;
  std::vector<dataPoint> m_data;
  Particles m_particles;
  EnergyType m_energyType = EnergyType::RIGIDITY;  // Track energy type for chi2 computation
};

class Chi2IH : public Chi2 {
 public:
  Chi2IH(const std::string& name, const Particles& particles, const double& phi) : Chi2(name, particles, phi) {
    constexpr double yunits = 1. / (CGS::GeV * CGS::m2 * CGS::sec * CGS::sr);
    readKISSfile("data/H_AMS-02_R.txt", xunits, yunits);
  }

 protected:
  double getModel(const double& R, const double& phi) const override;
  itParticle itH1 = findParticle(H1);
  itParticle itH2 = findParticle(H2);
  itParticle itH1_ter = findParticle(H1_ter);
};

class Chi2IHe : public Chi2 {
 public:
  Chi2IHe(const std::string& name, const Particles& particles, const double& phi) : Chi2(name, particles, phi) {
    constexpr double yunits = 1. / (CGS::GeV * CGS::m2 * CGS::sec * CGS::sr);
    readKISSfile("data/He_AMS-02_R.txt", xunits, yunits);
  }

 protected:
  double getModel(const double& R, const double& phi) const override;
  itParticle itHe3 = findParticle(He3);
  itParticle itHe4 = findParticle(He4);
};

class Chi2BC : public Chi2 {
 public:
  Chi2BC(const std::string& name, const Particles& particles, const double& phi) : Chi2(name, particles, phi) {
    readKISSfile("data/BC_AMS-02_R.txt", xunits, 1.);
  }

 protected:
  double getModel(const double& R, const double& phi) const override;
  itParticle itB10 = findParticle(B10);
  itParticle itB11 = findParticle(B11);
  itParticle itC12 = findParticle(C12);
  itParticle itC13 = findParticle(C13);
  itParticle itC14 = findParticle(C14);
};

// ============================================================================
// Generic Chi2Composite: Flexible model for combining isotopes with rules
// ============================================================================

// Structure to specify an isotope group
struct IsotopeGroup {
  std::vector<PID> pids;           // List of PIDs in this group
  double scaling = 1.0;            // Overall scaling factor for this group
  bool perNucleon = false;         // If true, divide by A (mass number)
};

// Generic composite model class
class Chi2Composite : public Chi2 {
 public:
  Chi2Composite(const std::string& name, 
                const Particles& particles,
                const double& phi,
                CompositionMode mode,
                EnergyType energyType,
                const std::string& dataFile,
                const double& xunits,
                const double& yunits,
                const std::vector<IsotopeGroup>& numeratorGroups,
                const std::vector<IsotopeGroup>& denominatorGroups = {})
      : Chi2(name, particles, phi),
        m_mode(mode),
        m_numeratorGroups(numeratorGroups),
        m_denominatorGroups(denominatorGroups),
        m_isManyPar(false) {
    // Store energy type in base class for computeChi2 to use
    this->m_energyType = energyType;
    
    // Load data from file - works for all energy types
    if (!dataFile.empty()) {
      readCRDBfile(dataFile, xunits, yunits);
    }
    
    // Cache particle iterators for all referenced PIDs
    for (const auto& group : numeratorGroups) {
      for (const auto& pid : group.pids) {
        m_cachedIterators[pid] = findParticle(pid);
      }
    }
    for (const auto& group : denominatorGroups) {
      for (const auto& pid : group.pids) {
        m_cachedIterators[pid] = findParticle(pid);
      }
    }
  }

  virtual ~Chi2Composite() = default;

 protected:
  double getModel(const double& energyVal, const double& phi) const override;
  
 private:
  // Helper to evaluate a group of isotopes
  double evaluateGroup(const std::vector<IsotopeGroup>& groups,
                       const double& energyVal,
                       const double& phi) const;

  CompositionMode m_mode;
  // NOTE: m_energyType is inherited from Chi2 base class - do NOT redeclare here
  std::vector<IsotopeGroup> m_numeratorGroups;
  std::vector<IsotopeGroup> m_denominatorGroups;
  bool m_isManyPar;  // For future: support complex energy dependencies
  
  // Cache for particle iterators
  mutable std::map<PID, itParticle> m_cachedIterators;
};

// ============================================================================
// Factory functions for creating Chi2Composite models
// ============================================================================

// Find all isotopes of a given element (Z) from the predefined PID list
std::vector<PID> findIsotopesOfElement(int Z);

// Get element symbol from atomic number (for auto data file finding)
std::string getElementSymbol(int Z);

// OVERLOADED: get_Chi2() - works with both manual and automatic isotope selection
// Manual version: pass vector of PID isotopes
std::unique_ptr<Chi2> get_Chi2(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    const std::vector<PID>& isotopes,
    const std::string& dataFile = "",
    const double& xunits = CGS::GeV,
    const double& yunits = 1.0 / (CGS::GeV * CGS::m2 * CGS::sec * CGS::sr));

// Automatic version: pass atomic number (int Z) to auto-detect all isotopes
std::unique_ptr<Chi2> get_Chi2(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    int Z,  // Atomic number - triggers automatic isotope detection
    const std::string& dataFile = "",
    const double& xunits = CGS::GeV,
    const double& yunits = 1.0 / (CGS::GeV * CGS::m2 * CGS::sec * CGS::sr));

// OVERLOADED: get_Chi2_ratio() - works with both manual and automatic isotope selection
// Manual version: pass vectors of PID isotopes for numerator and denominator
std::unique_ptr<Chi2> get_Chi2_ratio(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    const std::vector<PID>& numeratorIsotopes,
    const std::vector<PID>& denominatorIsotopes,
    const std::string& dataFile = "",
    const double& xunits = CGS::GeV,
    const double& yunits = 1.0);

// Automatic version: pass atomic numbers (int Z) for numerator and denominator
std::unique_ptr<Chi2> get_Chi2_ratio(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    int numeratorZ,     // Atomic number for numerator - triggers automatic detection
    int denominatorZ,   // Atomic number for denominator - triggers automatic detection
    const std::string& dataFile = "",
    const double& xunits = CGS::GeV,
    const double& yunits = 1.0);

// ============================================================================
// DAMPE-specific factory functions for flux and ratio data
// ============================================================================

// Create Chi2 for DAMPE total flux data (as function of total kinetic energy Ek in GeV)
std::unique_ptr<Chi2> get_Chi2_DAMPE_flux(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    int Z,  // Atomic number for element
    const std::string& dataFile = "",
    const double& xunits = CGS::GeV,
    const double& yunits = 1.0 / (CGS::GeV * CGS::m2 * CGS::sec * CGS::sr));

// Create Chi2 for DAMPE ratio data (as function of kinetic energy per nucleon Ekn in GeV)
std::unique_ptr<Chi2> get_Chi2_DAMPE_ratio(
    const std::string& name,
    const Particles& particles,
    const double& phi,
    int numeratorZ,
    int denominatorZ,
    const std::string& dataFile = "",
    const double& xunits = CGS::GeV,
    const double& yunits = 1.0);

}  // namespace CRAMS

// class Chi2_B: public Chi2 {
// public:
// 	Chi2_B(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		constexpr double units = 1. / (cgs::GeV * cgs::m2 * cgs::sec * cgs::sr);
// 		read_datafile("data/B_AMS02_rig.txt", units);
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_B10 = find_ptr(B10);
// 	ptr_Particle ptr_B11 = find_ptr(B11);
// };

// class Chi2_C: public Chi2 {
// public:
// 	Chi2_C(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		constexpr double units = 1. / (cgs::GeV * cgs::m2 * cgs::sec * cgs::sr);
// 		read_datafile("data/C_AMS02_rig.txt", units);
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_C12 = find_ptr(C12);
// 	ptr_Particle ptr_C13 = find_ptr(C13);
// 	ptr_Particle ptr_C14 = find_ptr(C14);
// };

// class Chi2_N: public Chi2 {
// public:
// 	Chi2_N(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		constexpr double units = 1. / (cgs::GeV * cgs::m2 * cgs::sec * cgs::sr);
// 		read_datafile("data/N_AMS02_rig.txt", units);
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_N14 = find_ptr(N14);
// 	ptr_Particle ptr_N15 = find_ptr(N15);
// };

// class Chi2_O: public Chi2 {
// public:
// 	Chi2_O(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		constexpr double units = 1. / (cgs::GeV * cgs::m2 * cgs::sec * cgs::sr);
// 		read_datafile("data/O_AMS02_rig.txt", units);
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_O16 = find_ptr(O16);
// 	ptr_Particle ptr_O17 = find_ptr(O17);
// 	ptr_Particle ptr_O18 = find_ptr(O18);
// };

// class Chi2_HeO: public Chi2 {
// public:
// 	Chi2_HeO(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		read_datafile("data/HeO_AMS02_rig.txt");
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_He3 = find_ptr(He3);
// 	ptr_Particle ptr_He4 = find_ptr(He4);
// 	ptr_Particle ptr_O16 = find_ptr(O16);
// 	ptr_Particle ptr_O17 = find_ptr(O17);
// 	ptr_Particle ptr_O18 = find_ptr(O18);
// };

// class Chi2_BeB: public Chi2 {
// public:
// 	Chi2_BeB(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		read_datafile("data/BeB_AMS02_rig.txt");
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_Be7 = find_ptr(Be7);
// 	ptr_Particle ptr_Be9 = find_ptr(Be9);
// 	ptr_Particle ptr_Be10 = find_ptr(Be10);
// 	ptr_Particle ptr_B10 = find_ptr(B10);
// 	ptr_Particle ptr_B11 = find_ptr(B11);
// };

// class Chi2_BeB_statsonly: public Chi2 {
// public:
// 	Chi2_BeB_statsonly(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		read_datafile_statsonly("data/BeB_AMS02_rig.txt");
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_Be7 = find_ptr(Be7);
// 	ptr_Particle ptr_Be9 = find_ptr(Be9);
// 	ptr_Particle ptr_Be10 = find_ptr(Be10);
// 	ptr_Particle ptr_B10 = find_ptr(B10);
// 	ptr_Particle ptr_B11 = find_ptr(B11);
// };

// class Chi2_BeC: public Chi2 {
// public:
// 	Chi2_BeC(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		read_datafile("data/BeC_AMS02_rig.txt");
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_Be7 = find_ptr(Be7);
// 	ptr_Particle ptr_Be9 = find_ptr(Be9);
// 	ptr_Particle ptr_Be10 = find_ptr(Be10);
// 	ptr_Particle ptr_C12 = find_ptr(C12);
// 	ptr_Particle ptr_C13 = find_ptr(C13);
// 	ptr_Particle ptr_C14 = find_ptr(C14);
// };

// class Chi2_BeO: public Chi2 {
// public:
// 	Chi2_BeO(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		read_datafile("data/BeO_AMS02_rig.txt");
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_Be7 = find_ptr(Be7);
// 	ptr_Particle ptr_Be9 = find_ptr(Be9);
// 	ptr_Particle ptr_Be10 = find_ptr(Be10);
// 	ptr_Particle ptr_O16 = find_ptr(O16);
// 	ptr_Particle ptr_O17 = find_ptr(O17);
// 	ptr_Particle ptr_O18 = find_ptr(O18);
// };

// class Chi2_BO: public Chi2 {
// public:
// 	Chi2_BO(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		read_datafile("data/BO_AMS02_rig.txt");
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_B10 = find_ptr(B10);
// 	ptr_Particle ptr_B11 = find_ptr(B11);
// 	ptr_Particle ptr_O16 = find_ptr(O16);
// 	ptr_Particle ptr_O17 = find_ptr(O17);
// 	ptr_Particle ptr_O18 = find_ptr(O18);
// };

// class Chi2_CO: public Chi2 {
// public:
// 	Chi2_CO(const Particles& particles, const double& phi) :
// 			Chi2(particles, phi) {
// 		read_datafile("data/CO_AMS02_rig.txt");
// 	}
// protected:
// 	double get_model(const double& R, const double& phi) const override;
// 	ptr_Particle ptr_C12 = find_ptr(C12);
// 	ptr_Particle ptr_C13 = find_ptr(C13);
// 	ptr_Particle ptr_C14 = find_ptr(C14);
// 	ptr_Particle ptr_O16 = find_ptr(O16);
// 	ptr_Particle ptr_O17 = find_ptr(O17);
// 	ptr_Particle ptr_O18 = find_ptr(O18);
// };

// ============================================================================
// Utility functions for data file discovery
// ============================================================================

#include <fstream>
#include <sys/stat.h>

namespace CRAMS {

// Helper to check if file exists
inline bool fileExists(const std::string& path) {
  struct stat buffer;
  return (stat(path.c_str(), &buffer) == 0);
}

// Search for a data file by element name, checking multiple naming conventions
// Looks for files in the data/ directory (relative to executable or absolute path)
inline std::string findDataFile(const std::string& elementName) {
  // Common naming patterns to check (in order of preference)
  std::vector<std::string> patterns = {
      elementName + "_AMS-02_R.txt",           // Short form: H_AMS-02_R.txt
      "AMS-02_" + elementName + "_rigidity.txt" // Long form: AMS-02_H_rigidity.txt
  };
  
  // Possible data directory locations
  std::vector<std::string> searchPaths = {
      "data",                                   // Relative to executable
      "./data",
      "../data",
      "../../data"
  };
  
  // Search for the file
  for (const auto& searchPath : searchPaths) {
    for (const auto& pattern : patterns) {
      std::string filepath = searchPath + "/" + pattern;
      if (fileExists(filepath)) {
        std::cout << "Found data file: " << filepath << std::endl;
        return filepath;
      }
    }
  }
  
  // File not found - return empty string (caller will skip data loading)
  return "";
}

}  // namespace CRAMS

#endif /* INCLUDE_CHI2_H_ */
