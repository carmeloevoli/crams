#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

#include "crams/core/cgs.h"
#include "crams/inelastic.h"
#include "crams/particle.h"
#include "crams/physics/grammage.h"
#include "crams/physics/losses.h"

namespace {

constexpr size_t kIntegrationLimit = 2000;
constexpr double kEpsRel = 1e-3;
constexpr double kEnergyIntegrationRatio = 1e3;
constexpr double kMinAttenuation = 1e-14;

using GslWorkspace = std::unique_ptr<gsl_integration_workspace, decltype(&gsl_integration_workspace_free)>;

double computeIntegralQags(gsl_integration_workspace* workspace, gsl_function* function, double xLow, double xHigh) {
  double result = 0.;
  double error = 0.;
  gsl_integration_qags(function, xLow, xHigh, 0., kEpsRel, kIntegrationLimit, workspace, &result, &error);
  return result;
}

double computeIntegralQag(gsl_integration_workspace* workspace, gsl_function* function, double xLow, double xHigh) {
  double result = 0.;
  double error = 0.;
  constexpr int integrationRule = 2;
  gsl_integration_qag(function, xLow, xHigh, 0., kEpsRel, kIntegrationLimit, integrationRule, workspace, &result,
                      &error);
  return result;
}

double geometricMidpoint(double left, double right) { return std::sqrt(left * right); }

struct TransportCoefficients {
  double attenuation = 0.;
  double source = 0.;
};

TransportCoefficients transportCoefficientsAt(const CRAMS::Particle& particle, double energy) {
  const double lambda1 = particle.Lambda_1(energy);
  const double lambda2 = particle.Lambda_2(energy);
  const double source = particle.Q_total(energy);

  if (!std::isfinite(lambda1) || !std::isfinite(lambda2) || !std::isfinite(source) || lambda2 <= 0.)
    throw std::runtime_error(
        "Particle: invalid transport coefficient in numerical solver: lambda1 = " + std::to_string(lambda1) +
        "; lambda2 = " + std::to_string(lambda2) + "; source = " + std::to_string(source));

  return {lambda1 / lambda2, source / lambda2};
}

std::vector<TransportCoefficients> buildTransportCoefficients(const CRAMS::Particle& particle,
                                                              const std::vector<double>& energies) {
  std::vector<TransportCoefficients> coefficients(energies.size());

  for (size_t i = 0; i < energies.size(); ++i) {
    coefficients[i] = transportCoefficientsAt(particle, energies[i]);
  }

  return coefficients;
}

double crankNicolsonStep(double nextIntensity, double energyStep, const TransportCoefficients& current,
                         const TransportCoefficients& next) {
  if (!std::isfinite(energyStep) || energyStep <= 0.)
    throw std::runtime_error("Particle: energy grid must be strictly increasing");

  const double halfStep = 0.5 * energyStep;
  const double numerator =
      (1. - halfStep * next.attenuation) * nextIntensity + halfStep * (current.source + next.source);
  const double denominator = 1. + halfStep * current.attenuation;

  if (!std::isfinite(numerator) || !std::isfinite(denominator) || denominator == 0.)
    throw std::runtime_error("Particle: singular Crank-Nicolson step");

  const double intensity = numerator / denominator;
  if (!std::isfinite(intensity)) throw std::runtime_error("Particle: non-finite Crank-Nicolson intensity");

  return intensity;
}

}  // namespace

namespace CRAMS {

double Particle::Lambda_1(double T) const {
  if (!m_X) throw std::runtime_error("Particle: buildGrammage must be called before computeIntensity");
  if (!m_dEdX) throw std::runtime_error("Particle: buildLosses must be called before computeIntensity");

  const double escape = 1. / m_X->get(T);
  const double inelastic = (m_sigmaIn) ? m_sigmaIn->getXsecOnISM(m_pid, T) / CGS::meanISMmass : 0.;
  const double lossesDerivative = m_dEdX->getDerivative(T);
  return escape + inelastic + lossesDerivative;
}

double Particle::Lambda_2(double T) const {
  if (!m_dEdX) throw std::runtime_error("Particle: buildLosses must be called before computeIntensity");
  return std::fabs(m_dEdX->get(T));
}

double Particle::internalIntegrand(double T_second) const { return Lambda_1(T_second) / Lambda_2(T_second); }

double gslParticleClassExpWrapper(double x, void* params) {
  const double TSecond = std::exp(x);
  const auto particle = static_cast<const Particle*>(params);
  return TSecond * particle->internalIntegrand(TSecond);
}

double Particle::ExpIntegral(double T, double T_prime) const {
  gsl_function F;
  F.params = const_cast<Particle*>(this);
  F.function = &gslParticleClassExpWrapper;

  GslWorkspace workspace(gsl_integration_workspace_alloc(kIntegrationLimit), gsl_integration_workspace_free);
  return computeIntegralQags(workspace.get(), &F, std::log(T), std::log(T_prime));
}

double Particle::externalIntegrand(double T_prime, double T) const {
  double value = Q_total(T_prime) * std::exp(-ExpIntegral(T, T_prime));
  value /= Lambda_2(T_prime);
  return value;
}

struct GslParticleParams {
  double T;
  const Particle* particle;
};

double gslParticleClassWrapper(double x, void* params) {
  const auto p = static_cast<const GslParticleParams*>(params);
  const double TPrime = std::exp(x);
  return TPrime * p->particle->externalIntegrand(TPrime, p->T);
}

double Particle::computeFluxAtEnergy(double T) const {
  GslParticleParams params = {T, this};
  gsl_function F;
  F.params = &params;
  F.function = &gslParticleClassWrapper;

  GslWorkspace workspace(gsl_integration_workspace_alloc(kIntegrationLimit), gsl_integration_workspace_free);
  return computeIntegralQag(workspace.get(), &F, std::log(T), std::log(kEnergyIntegrationRatio * T));
}

void Particle::computeFluxAtEnergyCrankNicolson() {
  if (m_T.size() < 2) return;

  const auto coefficients = buildTransportCoefficients(*this, m_T);
  const size_t last = m_T.size() - 1;
  m_I_T[last] = 0.;

  for (size_t i = last; i-- > 0;) {
    const double energyStep = m_T[i + 1] - m_T[i];
    m_I_T[i] = crankNicolsonStep(m_I_T[i + 1], energyStep, coefficients[i], coefficients[i + 1]);
  }
}

void Particle::computeFluxAtEnergyExponential() {
  if (m_T.size() < 2) return;

  const size_t last = m_T.size() - 1;
  m_I_T[last] = 0.;

  for (size_t i = last; i-- > 0;) {
    const double h = m_T[i + 1] - m_T[i];
    const double T_mid = geometricMidpoint(m_T[i], m_T[i + 1]);
    const double lambda1 = Lambda_1(T_mid);
    const double lambda2 = Lambda_2(T_mid);
    if (lambda2 <= 0.) throw std::runtime_error("Particle: Lambda_2 must be positive in numerical solver");

    const double source = Q_total(T_mid) / lambda2;
    const double attenuation = lambda1 / lambda2;

    if (std::fabs(attenuation) < kMinAttenuation) {
      m_I_T[i] = m_I_T[i + 1] + source * h;
      continue;
    }

    const double tau = attenuation * h;
    const double expFactor = std::exp(-tau);
    const double oneMinusExp = -std::expm1(-tau);
    m_I_T[i] = expFactor * m_I_T[i + 1] + source / attenuation * oneMinusExp;
  }
}

}  // namespace CRAMS
