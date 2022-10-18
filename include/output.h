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
  void dumpSpectraEkn() const;

 private:
  const Particles& m_particles;
  double m_phi = 0;
  size_t m_id = 0;
  std::vector<double> m_R;
  std::string m_simname;

 private:
  double getFluxChargeGroup(const int Z, const double& R) const;
  double getFluxChargeGroupEkn(const int Z, const double& T) const;
};

}  // namespace CRAMS

#endif  // INCLUDE_OUTPUT_MANAGER_H_
