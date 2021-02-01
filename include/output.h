#ifndef INCLUDE_OUTPUT_MANAGER_H_
#define INCLUDE_OUTPUT_MANAGER_H_

#include "input.h"
#include "particle.h"

namespace CRAMS {

class OutputManager {
 public:
  OutputManager(const Particles& particles, const Input& input);
  virtual ~OutputManager();
  void dumpSpectra() const;

 private:
  const Particles& m_particles;
  double m_phi = 0;
  size_t m_id = 0;
  std::vector<double> m_R;

 private:
  double getFluxChargeGroup(const int Z, const double& R) const;

  //   ptr_Particle ptr_H1_ter = find_ptr(H1_ter);
  //   ptr_Particle ptr_H1 = find_ptr(H1);
  //   ptr_Particle ptr_H2 = find_ptr(H2);
  //   ptr_Particle ptr_He3 = find_ptr(He3);
  //   ptr_Particle ptr_He4 = find_ptr(He4);
  //   ptr_Particle ptr_Li6 = find_ptr(Li6);
  //   ptr_Particle ptr_Li7 = find_ptr(Li7);
  //   ptr_Particle ptr_Be7 = find_ptr(Be7);
  //   ptr_Particle ptr_Be9 = find_ptr(Be9);
  //   ptr_Particle ptr_Be10 = find_ptr(Be10);
  //   ptr_Particle ptr_B10 = find_ptr(B10);
  //   ptr_Particle ptr_B11 = find_ptr(B11);
  //   ptr_Particle ptr_C12 = find_ptr(C12);
  //   ptr_Particle ptr_C13 = find_ptr(C13);
  //   ptr_Particle ptr_C14 = find_ptr(C14);
  //   ptr_Particle ptr_N14 = find_ptr(N14);
  //   ptr_Particle ptr_N15 = find_ptr(N15);
  //   ptr_Particle ptr_O16 = find_ptr(O16);
  //   ptr_Particle ptr_O17 = find_ptr(O17);
  //   ptr_Particle ptr_O18 = find_ptr(O18);
  //   ptr_Particle ptr_Fe56 = find_ptr(Fe56);
  //   ptr_Particle find_ptr(const PID& pid);

  //   std::string spectra_filename;
  //   std::string ratios_filename;

  //   double H(const double& R) const;
  //   double He(const double& R) const;
  //   double Li(const double& R) const;
  //   double Be(const double& R) const;
  //   double B(const double& R) const;
  //   double C(const double& R) const;
  //   double N(const double& R) const;
  //   double O(const double& R) const;
  //   double Fe(const double& R) const;

  //   double Be_ratio(const double& R) const;
  //   double He_ratio(const double& R) const;
};

}  // namespace CRAMS

#endif  // INCLUDE_OUTPUT_MANAGER_H_
