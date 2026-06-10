#include "crams/primary.h"

#include <plog/Log.h>

#include <cmath>

#include "crams/core/cgs.h"
#include "crams/utils/numeric.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

using Utilities::pow2;

PrimarySource::PrimarySource(const PID& pid, double abundance, double slope, double surfaceDensity)
    : m_pid(pid), m_slope(slope) {
  if (abundance > 0.) {
    static constexpr double L_SN_surface =
        CGS::E_SN * CGS::snRate / M_PI / pow2(CGS::galaxySize);
    m_norm = static_cast<double>(pid.getA()) * abundance * L_SN_surface
             / (surfaceDensity * Numeric::gammaIntegral(slope) * pow2(CGS::protonMassC2));
  }
}

PrimarySource::~PrimarySource() { LOGD << "deleted PrimarySource for particle " << m_pid; }

double PrimarySource::get(double T) const {
  if (m_norm <= 0.) return 0.;
  const double pc = Utilities::T2pc(T, m_pid);
  return m_norm / Utilities::T2beta(T) * std::pow(pc / CGS::protonMassC2, 2. - m_slope);
}

}  // namespace CRAMS
