#include "crams/physics/grammage.h"

#include <plog/Log.h>

#include <cmath>

#include "crams/core/cgs.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

using Utilities::pow2;

namespace {

const double kLn2 = std::log(2.);

double halfLifeToMeanLifetime(double halfLife) { return halfLife / kLn2; }

}  // namespace

Grammage::Grammage(const PID& pid, const Input& input) : m_pid(pid) { setParameters(input); }

Grammage::Grammage(const PID& pid, const Input& input, double decayHalfLifeAtRest)
    : m_pid(pid), m_decayHalfLifeAtRest(decayHalfLifeAtRest) {
  setParameters(input);
}

Grammage::~Grammage() { LOGD << "deleted Grammage for particle " << m_pid; }

void Grammage::setParameters(const Input& input) {
  m_norm = input.mu() * CGS::cLight / 2. / input.v_A();
  m_v_A = input.v_A();
  m_H = input.H();
  m_D_0 = input.D_0();
  m_R_b = input.R_b();
  m_delta = input.delta();
  m_ddelta = input.ddelta();
  m_smoothness = input.smoothness();
}

double Grammage::D(double T) const {
  const double R = Utilities::T2pc(T, m_pid) / std::abs(m_pid.getZ());
  const double x = R / m_R_b;
  const double smooth = std::pow(1. + std::pow(x, m_ddelta / m_smoothness), m_smoothness);
  return m_D_0 * Utilities::T2beta(T) * std::pow(R / CGS::GeV, m_delta) / smooth + 2. * m_v_A * m_H;
}

double Grammage::get(double T) const {
  const double beta = Utilities::T2beta(T);
  const double d = D(T);
  const double escape_depth = m_v_A * m_H / d;

  if (m_decayHalfLifeAtRest < 0.) {
    // Stable particle: standard slab grammage
    return beta * m_norm * (1. - std::exp(-escape_depth));
  } else {
    // Unstable particle: decay modifies escape probability
    const double meanLifetimeAtRest = halfLifeToMeanLifetime(m_decayHalfLifeAtRest);
    const double tau_d = Utilities::T2gamma(T) * meanLifetimeAtRest;
    const double Delta = std::sqrt(1. + 4. * d / (pow2(m_v_A) * tau_d));
    const double exp_term = std::exp(-escape_depth * Delta);
    return beta * m_norm * 2. * (1. - exp_term) / ((1. + Delta) - (1. - Delta) * exp_term);
  }
}

}  // namespace CRAMS
