#ifndef CRAMS_SECONDARY_H_
#define CRAMS_SECONDARY_H_

#include <vector>

#include "crams/core/pid.h"

namespace CRAMS {

class SecondarySource {
 public:
  SecondarySource(const PID& pid, const std::vector<double>& T, const std::vector<double>& Q);
  ~SecondarySource();

  // Returns the interpolated secondary source at kinetic energy T, or 0 outside the grid.
  double get(double T) const;

 private:
  PID m_pid;
  std::vector<double> m_T;
  std::vector<double> m_Q;
};

}  // namespace CRAMS

#endif  // CRAMS_SECONDARY_H_
