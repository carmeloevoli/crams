#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include "particle.h"

namespace CRAMS {

#define LIMIT 2000
#define EPSREL 1e-3
#define INFTY4INTEGRAL 1e3

double computeIntegralQags(gsl_integration_workspace* w, gsl_function* F, double x_lo, double x_hi) {
  double result, error;
  gsl_integration_qags(F, x_lo, x_hi, 0, EPSREL, LIMIT, w, &result, &error);
  return result;
}

double computeIntegralQag(gsl_integration_workspace* w, gsl_function* F, double x_lo, double x_hi) {
  double result, error;
  int key = 2;
  gsl_integration_qag(F, x_lo, x_hi, 0, EPSREL, LIMIT, key, w, &result, &error);
  return result;
}

double Particle::Lambda_1(const double& T) const {
  return 1. / m_X->get(T) + m_sigmaIn->getXsecOnISM(T) / CGS::meanISMmass + m_dEdX->getDerivative(T);
}

double Particle::Lambda_2(const double& T) const {
  double value = m_dEdX->get(T);
  return std::fabs(value);
}

double Particle::internalIntegrand(const double& T_second) {
  double value = Lambda_1(T_second) / Lambda_2(T_second);
  return value;
}

double gslParticleClassExpWrapper(double x, void* pp) {
  double E_second = std::exp(x);
  Particle* particle = (Particle*)pp;
  return E_second * particle->internalIntegrand(E_second);
}

double Particle::ExpIntegral(const double& T, const double& T_prime) {
  double result = 0;
  gsl_function F;
  F.params = this;
  F.function = &gslParticleClassExpWrapper;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(LIMIT);
  result = computeIntegralQags(w, &F, std::log(T), std::log(T_prime));
  gsl_integration_workspace_free(w);
  return result;
}

double Particle::externalIntegrand(const double& T_prime, const double& T) {
  double value = Q_total(T_prime) * std::exp(-ExpIntegral(T, T_prime));
  value /= Lambda_2(T_prime);
  return value;
}

struct gsl_f_pars {
  double T;
  Particle* pt_Particle;
};

double gslParticleClassWrapper(double x, void* pp) {
  gsl_f_pars* p = (gsl_f_pars*)pp;
  double T = p->T;
  double T_prime = std::exp(x);
  return T_prime * p->pt_Particle->externalIntegrand(T_prime, T);
}

double Particle::computeFluxAtEnergy(const double& T) {
  double result = 0;
  gsl_f_pars pars = {T, this};
  gsl_function F;
  F.params = &pars;
  F.function = &gslParticleClassWrapper;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(LIMIT);
  result = computeIntegralQag(w, &F, std::log(T), std::log(INFTY4INTEGRAL * T));
  gsl_integration_workspace_free(w);
  return result;
}

}  // namespace CRAMS