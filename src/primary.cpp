#include "primary.h"

#include <plog/Log.h>

#include "cgs.h"
#include "utilities.h"

namespace CRAMS {

PrimarySource::PrimarySource() {}

PrimarySource::~PrimarySource() { LOGD << "deleted PrimarySource for particle " << m_pid; }

PrimarySource::PrimarySource(const PID& pid, const double& abundance, const double& slope, const double& surfaceDensity)
    : m_pid(pid), m_slope(slope) {
  if (abundance > 0.) {
    constexpr double L_SN_surface = CGS::E_SN * CGS::snRate / M_PI / pow2(CGS::galaxySize);
    m_energyFactor = (double)pid.getA() * abundance * L_SN_surface;
    m_energyFactor /= surfaceDensity * Utilities::GammaIntegral(slope) * pow2(CGS::protonMassC2);
  }
}

double PrimarySource::get(const double& T) const {
  double value = 0.;
  if (m_energyFactor > 0.) {
    const double beta = Utilities::T2beta(T);
    const double pc = Utilities::T2pc(T, m_pid);
    value = m_energyFactor / beta;
    value *= std::pow(pc / CGS::protonMassC2, 2. - m_slope);
  }
  return value;
}

}  // namespace CRAMS