#include "crams/physics/losses.h"

#include <gsl/gsl_deriv.h>
#include <plog/Log.h>

#include <cmath>

#include "crams/core/cgs.h"
#include "crams/utils/utilities.h"

namespace {

double lossesWrapper(double x, void* pp) {
  const auto* losses = static_cast<const CRAMS::Losses*>(pp);
  return losses->dEdX_adiabatic(x) + losses->dEdX_ionization(x);
}

}  // namespace

namespace CRAMS {

using Utilities::pow2;

Losses::Losses(const PID& pid, const Input& input) : m_pid(pid) {
  m_factorAdv = 2. * input.v_A() / 3. / input.mu() / CGS::cLight;
}

Losses::~Losses() { LOGD << "deleted Losses for particle " << m_pid; }

double Losses::betheBlochLog(double T) const {
  const double beta = Utilities::T2beta(T);
  const double gamma = Utilities::T2gamma(T);
  const double mA = static_cast<double>(m_pid.getA()) * CGS::protonMass;
  const double beta2 = pow2(beta);
  const double gamma2 = pow2(gamma);
  const double Q_max = 2. * CGS::electronMassC2 * beta2 * gamma2 / (1. + 2. * gamma * CGS::electronMass / mA);
  const double arg = 2. * CGS::electronMassC2 * (gamma2 - 1.) * Q_max;
  const double B_H = std::log(arg / pow2(CGS::IsH)) - 2. * beta2;
  const double B_He = std::log(arg / pow2(CGS::IsHe)) - 2. * beta2;
  return B_H + CGS::f_He * B_He;
}

double Losses::dEdX_adiabatic(double T) const { return -m_factorAdv * std::sqrt(T * (T + 2. * CGS::protonMassC2)); }

double Losses::dEdX_ionization(double T) const {
  constexpr double k = 2. * M_PI * pow2(CGS::electronRadius) * CGS::electronMassC2;
  const double Z = static_cast<double>(m_pid.getZ());
  const double A = static_cast<double>(m_pid.getA());
  const double beta2 = pow2(Utilities::T2beta(T));
  return -k * Z * Z * betheBlochLog(T) / (CGS::protonMass * (1. + 4. * CGS::f_He) * A * beta2);
}

double Losses::dTdt_ionization(double T, double n_H) const {
  constexpr double k = 2. * M_PI * pow2(CGS::electronRadius) * CGS::cLight * CGS::electronMassC2;
  const double Z = static_cast<double>(m_pid.getZ());
  const double A = static_cast<double>(m_pid.getA());
  return k * Z * Z / A * Utilities::T2beta(T) * n_H * betheBlochLog(T);
}

double Losses::get(double T) const { return dEdX_adiabatic(T) + dEdX_ionization(T); }

double Losses::getDerivative(double T) const {
  gsl_function F;
  F.function = &lossesWrapper;
  F.params = const_cast<Losses*>(this);  // GSL needs void*; only const methods are called
  double result, abserr;
  gsl_deriv_central(&F, T, 0.01 * T, &result, &abserr);
  return result;
}

}  // namespace CRAMS
