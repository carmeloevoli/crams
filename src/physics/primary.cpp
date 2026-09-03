#include "crams/physics/primary.h"

#include <plog/Log.h>

#include <cmath>

#include "crams/core/cgs.h"
#include "crams/utils/numeric.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

using Utilities::pow2;

double SpectralBreak::computeModifier(double T) const {
  return std::pow(1. + std::pow(T / m_T, 1. / m_omega), -(m_deltaSlope * m_omega));
}

double ErfcCutoff::computeModifier(double T) const {
  const double erfArg = (std::log10(T) - m_lgT - pow2(m_sigma) * m_beta_ln10) / (M_SQRT2 * m_sigma);
  if (m_isLower) {
    return 0.5 * (1 + std::erf(erfArg));
  } else {
    return 0.5 * std::erfc(erfArg);
  }
}

double ExpCutoff::computeModifier(double T) const {
  if (m_isLower) {
    return std::exp(-(m_T / T) * m_Delta);
  } else {
    return std::exp(-(T / m_T) * m_Delta);
  }
}

PrimarySource::PrimarySource(const PID& pid, double abundance, double slope, double surfaceDensity)
    : m_pid(pid), m_slope(slope) {
  if (abundance > 0.) {
    static constexpr double L_SN_surface = CGS::E_SN * CGS::snRate / M_PI / pow2(CGS::galaxySize);
    m_norm = static_cast<double>(pid.getA()) * abundance * L_SN_surface /
             (surfaceDensity * Numeric::gammaIntegral(slope) * pow2(CGS::protonMassC2));
  }
}

PrimarySource::~PrimarySource() { LOGD << "deleted PrimarySource for particle " << m_pid; }

double PrimarySource::get(double T) const {
  if (m_norm <= 0.) return 0.;
  const double pc = Utilities::T2pc(T, m_pid);

  double spectrumModifier = 1.0;
  for (const auto& feature : m_features) {
    spectrumModifier *= feature->computeModifier(T);
  }

  return m_norm / Utilities::T2beta(T) * std::pow(pc / CGS::protonMassC2, 2. - m_slope) * spectrumModifier;
}

}  // namespace CRAMS
