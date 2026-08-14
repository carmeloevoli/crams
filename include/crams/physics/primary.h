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

  void resetSpectralFeatures();
  void setSpectralBreak(double R, double deltaSlope, double omega);
  void setErfCutoff(double R, double sigma, double beta);

 private:
  PID m_pid;
  double m_slope = 4.0;
  double m_norm = 0.;  // amplitude: SNR budget / (gas density × spectral integral) [/erg²]

  // source spectrum features
  double m_featureT = -1;  // <0 = featureless PL
  double m_featureLgT = -1;  // <0 = featureless PL

  // feature = smooth break
  bool m_isBreak = false;
  double m_breakDeltaSlope;
  double m_breakOmega;

  // feature = lognormal distribution of maximum energies of CR accelerators -> erf-shaped cutoff
  bool m_isLognorm = false;
  double m_lognormSigma;  // in decades
  double m_lognormBeta;   // PL index in the dependence of CR accelerator luminocity on the maximum energy.
                          // when convolving the individual cut-offs with population weight, we have
                          // W(Emax) \propto Emax^beta, beta~1 for standard models of SNR acceleration
};

}  // namespace CRAMS

#endif  // CRAMS_PHYSICS_PRIMARY_H_
