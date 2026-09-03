#ifndef CRAMS_PHYSICS_PRIMARY_H_
#define CRAMS_PHYSICS_PRIMARY_H_

#include <memory>
#include <vector>

#include "crams/core/pid.h"

namespace CRAMS {

struct SpectralFeature {
  virtual ~SpectralFeature() = default;
  virtual double computeModifier(double T) const = 0;
};

struct SpectralBreak : public SpectralFeature {
  SpectralBreak(double T, double deltaSlope, double omega) : m_T(T), m_deltaSlope(deltaSlope), m_omega(omega) {};
  double computeModifier(double T) const override;

 private:
  double m_T;
  double m_deltaSlope;
  double m_omega;
};

struct ErfcCutoff : public SpectralFeature {
  ErfcCutoff(double T, double sigma, double beta, bool isLower)
      : m_lgT(std::log10(T)), m_sigma(sigma), m_beta_ln10(beta * std::log(10.)), m_isLower(isLower) {};
  double computeModifier(double T) const override;

 private:
  double m_lgT;
  double m_sigma;      // in decades
  double m_beta_ln10;  // PL index in the dependence of CR accelerator luminocity on the maximum energy.
                       // when convolving the individual cut-offs with population weight, we have
                       // W(Emax) \propto Emax^beta, beta~1 for standard models of SNR acceleration
                       // here we store beta * ln(10) to avoid recomputing it
  bool m_isLower;
};

struct ExpCutoff : public SpectralFeature {
  ExpCutoff(double T, double Delta, bool isLower) : m_T(T), m_Delta(Delta), m_isLower(isLower) {}
  double computeModifier(double T) const override;

 private:
  double m_T;
  double m_Delta;
  bool m_isLower;
};

class PrimarySource {
 public:
  PrimarySource(const PID& pid, double abundance, double slope, double surfaceDensity);
  PrimarySource(const PrimarySource&) = delete;
  PrimarySource& operator=(const PrimarySource&) = delete;
  PrimarySource(PrimarySource&&) noexcept = default;
  PrimarySource& operator=(PrimarySource&&) noexcept = default;
  ~PrimarySource();

  // Differential primary source at kinetic energy per nucleon T
  double get(double T) const;

  void addFeature(std::shared_ptr<SpectralFeature> feature) { m_features.push_back(std::move(feature)); }

 private:
  PID m_pid;
  double m_slope = 4.0;
  double m_norm = 0.;  // amplitude: SNR budget / (gas density × spectral integral) [/erg²]

  // multiplicative source spectrum features
  std::vector<std::shared_ptr<SpectralFeature>> m_features;
};

}  // namespace CRAMS

#endif  // CRAMS_PHYSICS_PRIMARY_H_
