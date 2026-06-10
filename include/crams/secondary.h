#ifndef CRAMS_SECONDARY_H_
#define CRAMS_SECONDARY_H_

#include <vector>

#include "crams/core/pid.h"
#include "crams/utils/utilities.h"

namespace CRAMS {
class SecondarySource {
 public:
  SecondarySource();
  SecondarySource(const PID& pid, const std::vector<double>& T, const std::vector<double>& Q);
  virtual ~SecondarySource();
  double get(const double& T) const;

 private:
  std::vector<double> m_T;
  std::vector<double> m_Q;
  PID m_pid;
};

}  // namespace CRAMS

#endif  // CRAMS_SECONDARY_H_
