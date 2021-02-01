#include "grammage.h"

#include <plog/Log.h>

#include <cmath>

#include "utilities.h"

namespace CRAMS {

Grammage::Grammage() {}  // TODO add cout for building classes with logger

Grammage::Grammage(const PID& pid, const Input& input) : m_pid(pid) { setParameters(input); }

Grammage::Grammage(const PID& pid, const Input& input, const double& tauDecayAtRest)
    : m_pid(pid), m_tauDecayAtRest(tauDecayAtRest) {
  setParameters(input);
}

Grammage::~Grammage() { LOGD << "deleted Grammage for particle " << m_pid; }

void Grammage::setParameters(const Input& input) {
  m_factor = input.mu * CGS::cLight / 2. / input.v_A;
  m_v_A = input.v_A;
  m_H = input.H;
  m_D_0 = input.D_0;
  m_R_b = input.R_b;
  m_delta = input.delta;
  m_ddelta = input.ddelta;
  m_s = input.smoothness;
}

double Grammage::D(const double& T) const {
  const double pc = Utilities::T2pc(T, m_pid);
  const double R = pc / m_pid.getZ();
  const double x = R / m_R_b;
  double value = Utilities::T2beta(T) * std::pow(R / CGS::GeV, m_delta);
  value /= std::pow(1. + std::pow(x, m_ddelta / m_s), m_s);
  return m_D_0 * value + 2.0 * m_v_A * m_H;
}

double Grammage::get(const double& T) const {
  const double beta = Utilities::T2beta(T);
  if (m_tauDecayAtRest < 0.) {
    double value = beta * m_factor;
    value *= 1. - std::exp(-m_v_A * m_H / D(T));
    return value;
  } else {
    double tau_d = Utilities::T2gamma(T) * m_tauDecayAtRest;
    double xi = m_v_A * m_H / D(T);
    double Delta = std::sqrt(1. + 4. * D(T) / pow2(m_v_A) / tau_d);
    double value = beta * m_factor;
    value *= 2. * (1. - std::exp(-xi * Delta));
    value /= (1. + Delta) - (1. - Delta) * std::exp(-xi * Delta);
    return value;
  }
}

}  // namespace CRAMS