#include "losses.h"

#include <gsl/gsl_deriv.h>
#include <plog/Log.h>

#include <cmath>
#include <string>

#include "cgs.h"
#include "utilities.h"

namespace CRAMS {

Losses::Losses() {}

Losses::Losses(const PID& pid, const Input& input) : m_pid(pid), m_mu(input.mu) {
  m_factorAdv = 2. * input.v_A / 3. / input.mu / CGS::cLight;
}

Losses::~Losses() { LOGD << "deleted Losses for particle " << m_pid; }

double Losses::get(const double& T) const { return dEdX_adiabatic(T) + dEdX_ionization(T); }

double Losses::dEdX_adiabatic(const double& T) const {
  double value = m_factorAdv * std::sqrt(T * (T + 2. * CGS::protonMassC2));
  return -value;
}

double Losses::dEdX_ionization(const double& T) const {
  const double beta = Utilities::T2beta(T);
  const double gamma = Utilities::T2gamma(T);
  constexpr double mec2 = CGS::electronMassC2;
  constexpr double twoPiRe2 = 2. * M_PI * pow2(CGS::electronRadius);
  double Q_max = 2. * CGS::electronMassC2 * pow2(beta) * pow2(gamma);
  Q_max /= 1. + 2. * gamma * CGS::electronMass / ((double)m_pid.getA() * CGS::protonMass);
  double B_H = std::log(2. * mec2 * (pow2(gamma) - 1.) * Q_max / pow2(CGS::IsH));
  B_H -= 2. * pow2(beta);
  double B_He = std::log(2. * mec2 * (pow2(gamma) - 1.) * Q_max / pow2(CGS::IsHe));
  B_He -= 2. * pow2(beta);
  double dEdx = twoPiRe2 * mec2 * (double)pow2(m_pid.getZ()) * (B_H + B_He * CGS::f_He);
  dEdx /= CGS::protonMass * (1. + 4. * CGS::f_He) * (double)m_pid.getA() * pow2(beta);
  return -dEdx;
}

double Losses::dTdt_ionization(const double& T, const double& n_H) const {  // TODO COMBINE WITH dEdX!
  const double beta = Utilities::T2beta(T);
  const double gamma = Utilities::T2gamma(T);
  constexpr double mec2 = CGS::electronMassC2;
  constexpr double twoPiRe2Cmec2 = 2. * M_PI * pow2(CGS::electronRadius) * CGS::cLight * mec2;
  const double Zsquared = (double)pow2(m_pid.getZ());
  const double A = (double)m_pid.getA();
  double Q_max = 2. * CGS::electronMassC2 * pow2(beta) * pow2(gamma);
  Q_max /= 1. + 2. * gamma * CGS::electronMass / ((double)m_pid.getA() * CGS::protonMass);
  double B_H = std::log(2. * mec2 * (pow2(gamma) - 1.) * Q_max / pow2(CGS::IsH));
  B_H -= 2. * pow2(beta);
  double B_He = std::log(2. * mec2 * (pow2(gamma) - 1.) * Q_max / pow2(CGS::IsHe));
  B_He -= 2. * pow2(beta);
  return twoPiRe2Cmec2 * Zsquared / A * beta * n_H * (B_H + B_He * CGS::f_He);
}

double gslLossesClassWrapper(double x, void* pp) {
  Losses* losses = (Losses*)pp;
  return losses->dEdX_adiabatic(x) + losses->dEdX_ionization(x);
}

double Losses::getDerivative(const double& T) {
  double result, abserr;
  gsl_function F;
  F.function = &gslLossesClassWrapper;
  F.params = this;
  double h_start = 0.01 * T;
  gsl_deriv_central(&F, T, h_start, &result, &abserr);
  return result;
}

}  // namespace CRAMS