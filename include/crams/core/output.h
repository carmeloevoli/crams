#ifndef CRAMS_CORE_OUTPUT_H_
#define CRAMS_CORE_OUTPUT_H_

#include <string>
#include <vector>

#include "crams/core/input.h"
#include "crams/particle.h"

namespace CRAMS {

class OutputManager {
 public:
  OutputManager(const Particles& particles, const Input& input);
  ~OutputManager() = default;

  void dumpSpectraRigidity() const;
  void dumpSpectraEkn() const;
  void dumpIsotopes() const;

 private:
  double getFluxChargeGroup(int Z, double R) const;
  double getFluxChargeIsotope(int Z, int A, double R) const;
  double getFluxChargeGroupEkn(int Z, double T) const;

  const Particles& m_particles;
  double m_phi = 0;
  size_t m_id = 0;
  std::vector<double> m_R;
  std::string m_simname;
};

}  // namespace CRAMS

#endif  // CRAMS_CORE_OUTPUT_H_
