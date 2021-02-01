#ifndef INCLUDE_PRIMARY_H_
#define INCLUDE_PRIMARY_H_

#include "pid.h"

namespace CRAMS {
class PrimarySource {
 public:
  PrimarySource();
  PrimarySource(const PID& pid, const double& abundance, const double& slope, const double& surfaceDensity);
  virtual ~PrimarySource();
  double get(const double& T) const;

 protected:
  PID m_pid;
  double m_slope = 4.0;
  double m_energyFactor = 0.;
};

}  // namespace CRAMS

#endif  // INCLUDE_PRIMARY_H_
