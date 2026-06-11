#ifndef CRAMS_PHYSICS_GRAMMAGE_H_
#define CRAMS_PHYSICS_GRAMMAGE_H_

#include "crams/core/input.h"
#include "crams/core/pid.h"

namespace CRAMS {

class Grammage {
 public:
  Grammage(const PID& pid, const Input& input);
  Grammage(const PID& pid, const Input& input, double decayHalfLifeAtRest);
  ~Grammage();

  // Diffusion coefficient at kinetic energy T [cm²/s]
  double D(double T) const;

  // Mean grammage traversed at kinetic energy T [g/cm²]
  double get(double T) const;

  double diffusionTimescale(double T) const { return m_H * m_H / D(T); }
  double advectionTimescale() const { return m_H / m_v_A; }

 private:
  void setParameters(const Input& input);

  PID m_pid;
  double m_norm = 0;  // mu * c / (2 * v_A) [g/cm²]
  double m_v_A = 0;
  double m_H = 0;
  double m_D_0 = 0;
  double m_R_b = 0;
  double m_delta = 0;
  double m_ddelta = 0;
  double m_smoothness = 0;
  double m_decayHalfLifeAtRest = -1;  // negative sentinel → stable particle
};

}  // namespace CRAMS

#endif  // CRAMS_PHYSICS_GRAMMAGE_H_
