#ifndef CRAMS_PRIMARY_H_
#define CRAMS_PRIMARY_H_

#include "crams/core/pid.h"

namespace CRAMS {

class PrimarySource {
 public:
  PrimarySource(const PID& pid, double abundance, double slope, double surfaceDensity);
  ~PrimarySource();

  // Differential primary source at kinetic energy per nucleon T
  double get(double T) const;

 private:
  PID m_pid;
  double m_slope = 4.0;
  double m_norm = 0.;  // amplitude: SNR budget / (gas density × spectral integral) [/erg²]
};

}  // namespace CRAMS

#endif  // CRAMS_PRIMARY_H_
