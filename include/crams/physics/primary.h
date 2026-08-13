#ifndef CRAMS_PHYSICS_PRIMARY_H_
#define CRAMS_PHYSICS_PRIMARY_H_

#include "crams/core/pid.h"

namespace CRAMS {

class PrimarySource {
 public:
  PrimarySource(const PID& pid, double abundance, double slope, double surfaceDensity);
  ~PrimarySource();

  // Differential primary source at kinetic energy per nucleon T
  double get(double T) const;

  void setSpectralBreak(double R, double deltaSlope, double omega);

 private:
  PID m_pid;
  double m_slope = 4.0;
  double m_norm = 0.;  // amplitude: SNR budget / (gas density × spectral integral) [/erg²]

  // source spectrum features
  double m_featureT = -1;  // <0 = featureless PL

  // feature = smooth break
  double m_breakDeltaSlope = 0.0;
  double m_breakOmega = 0.05;
};

}  // namespace CRAMS

#endif  // CRAMS_PHYSICS_PRIMARY_H_
