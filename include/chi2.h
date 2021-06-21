#ifndef INCLUDE_CHI2_H_
#define INCLUDE_CHI2_H_

#include <algorithm>
#include <string>
#include <vector>

#include "particle.h"
#include "pid.h"

namespace CRAMS {

struct dataPoint {
  double R;
  double I;
  std::pair<double, double> I_err;
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

  itParticle findParticle(const PID& pid) {
    auto it = find(m_particles.begin(), m_particles.end(), Particle(pid));
    bool isPresent = !(it == m_particles.end());
    return itParticle(isPresent, it);
  }

 protected:
  virtual double getModel(const double& R, const double& phi) const { return 0; }
  void readKISSfile(const std::string& filename, const double& xunits, const double& yunits = 1.0);

 protected:
  const double xunits = CGS::GeV;
  std::string m_name;
  double m_chi2 = 0;
  double m_phi = 0;
  std::vector<dataPoint> m_data;
  Particles m_particles;
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

#endif /* INCLUDE_CHI2_H_ */
