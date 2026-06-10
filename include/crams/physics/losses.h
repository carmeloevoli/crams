#ifndef CRAMS_PHYSICS_LOSSES_H_
#define CRAMS_PHYSICS_LOSSES_H_

#include "crams/core/input.h"
#include "crams/core/pid.h"

namespace CRAMS {

class Losses {
 public:
  Losses(const PID& pid, const Input& input);
  ~Losses();

  double get(double T) const;
  double dEdX_adiabatic(double T) const;
  double dEdX_ionization(double T) const;
  double getDerivative(double T) const;

  // Ionization cooling rate [erg/s] at kinetic energy T and hydrogen number density n_H
  double dTdt_ionization(double T, double n_H) const;

 private:
  // Effective Bethe-Bloch logarithm: B_H + f_He * B_He
  double betheBlochLog(double T) const;

  PID m_pid;
  double m_factorAdv = 0;  // 2*v_A / (3*mu*c) [cm²/g]
};

}  // namespace CRAMS

#endif  // CRAMS_PHYSICS_LOSSES_H_
