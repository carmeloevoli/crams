#ifndef INCLUDE_GRAMMAGE_H_
#define INCLUDE_GRAMMAGE_H_

#include "input.h"
#include "pid.h"

namespace CRAMS {
// namespace CRAMS
class Grammage {
 public:
  Grammage();
  Grammage(const PID& pid, const Input& input);
  Grammage(const PID& pid, const Input& input, const double& tauDecayAtRest);
  virtual ~Grammage();

  double D(const double& T) const;
  double get(const double& T) const;
  double diffusionTimescale(const double& T) const { return m_H * m_H / D(T); }
  double advectionTimescale() const { return m_H / m_v_A; }

 protected:
  void setParameters(const Input& input);

 protected:
  PID m_pid;
  double m_factor = 0;
  double m_v_A = 0;
  double m_H = 0;
  double m_D_0 = 0;
  double m_R_b = 0;
  double m_delta = 0;
  double m_ddelta = 0;
  double m_s = 0;
  double m_tauDecayAtRest = -1;
};

}  // namespace CRAMS

#endif  // INCLUDE_GRAMMAGE_H_
