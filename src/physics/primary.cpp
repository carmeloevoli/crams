#include "crams/physics/primary.h"

#include <plog/Log.h>

#include <cmath>
#include <stdexcept>

#include "crams/core/cgs.h"
#include "crams/utils/numeric.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

using Utilities::pow2;

PrimarySource::PrimarySource(const PID& pid, double abundance, double slope, double surfaceDensity)
    : m_pid(pid), m_slope(slope) {
  if (abundance > 0.) {
    static constexpr double L_SN_surface = CGS::E_SN * CGS::snRate / M_PI / pow2(CGS::galaxySize);
    m_norm = static_cast<double>(pid.getA()) * abundance * L_SN_surface /
             (surfaceDensity * Numeric::gammaIntegral(slope) * pow2(CGS::protonMassC2));
  }
}

PrimarySource::~PrimarySource() { LOGD << "deleted PrimarySource for particle " << m_pid; }

void PrimarySource::resetSpectralFeatures() {
  m_featureT = -1.0;
  m_isBreak = false;
  m_isLognorm = false;
}

void PrimarySource::setSpectralBreak(double R, double deltaSlope, double omega) {
  resetSpectralFeatures();
  m_isBreak = true;
  m_featureT = Utilities::R2T(R, m_pid);
  m_breakDeltaSlope = deltaSlope;
  m_breakOmega = omega;
};

void PrimarySource::setErfCutoff(double R, double sigma, double beta) {
  resetSpectralFeatures();
  m_isLognorm = true;
  m_featureT = Utilities::R2T(R, m_pid);
  m_featureLgT = std::log10(m_featureT);
  m_lognormSigma = sigma;
  m_lognormBeta = beta;
};

double PrimarySource::get(double T) const {
  if (m_norm <= 0.) return 0.;
  const double pc = Utilities::T2pc(T, m_pid);

  double spectrumModifier = 1.0;
  if (m_featureT >= 0) {
    if (m_isBreak) {
      // applying break in the source spectrum
      spectrumModifier = std::pow(1 + std::pow(T / m_featureT, 1 / m_breakOmega), -(m_breakDeltaSlope * m_breakOmega));
    } else if (m_isLognorm) {
      spectrumModifier = 0.5 * std::erfc((std::log10(T) - m_featureLgT - pow2(m_lognormSigma) * m_lognormBeta) /
                                         (M_SQRT2 * m_lognormSigma));
    } else {
      throw std::runtime_error("feature T is set, but no particular feature is selected with a flag");
    }
  }

  return m_norm / Utilities::T2beta(T) * std::pow(pc / CGS::protonMassC2, 2. - m_slope) * spectrumModifier;
}

}  // namespace CRAMS
