#ifndef CRAMS_LOSSES_H_
#define CRAMS_LOSSES_H_

#include "crams/core/input.h"
#include "crams/core/pid.h"

namespace CRAMS {
class Losses {
 public:
  Losses();
  Losses(const PID& pid, const Input& input);
  virtual ~Losses();
  double get(const double& T) const;
  double dEdX_adiabatic(const double& T) const;
  double dEdX_ionization(const double& T) const;
  double getDerivative(const double& T);
  double dTdt_ionization(const double& T, const double& n_H) const;

 protected:
  PID m_pid;
  double m_factorAdv = 0;
  double m_mu = 0;
};

}  // namespace CRAMS

#endif  // CRAMS_LOSSES_H_
