#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include <cmath>
#include <cstddef>
#include <memory>

#include "crams/core/cgs.h"
#include "crams/inelastic.h"
#include "crams/particle.h"
#include "crams/physics/grammage.h"
#include "crams/physics/losses.h"

namespace {

constexpr size_t kIntegrationLimit = 2000;
constexpr double kEpsRel = 1e-3;
constexpr double kEnergyIntegrationRatio = 1e3;

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

}  // namespace

namespace CRAMS {

double Particle::Lambda_1(double T) const {
  return 1. / m_X->get(T) + m_sigmaIn->getXsecOnISM(T) / CGS::meanISMmass + m_dEdX->getDerivative(T);
}

double Particle::Lambda_2(double T) const {
  return std::fabs(m_dEdX->get(T));
}

double Particle::internalIntegrand(double T_second) const {
  return Lambda_1(T_second) / Lambda_2(T_second);
}

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

}  // namespace CRAMS
