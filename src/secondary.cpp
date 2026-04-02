#include "secondary.h"

#include <plog/Log.h>

#include "gsl.h"

namespace CRAMS {

SecondarySource::SecondarySource() {}  // TODO have spallation inside, too stupid as it is?

SecondarySource::SecondarySource(const PID& pid, const std::vector<double>& T, const std::vector<double>& Q)
    : m_pid(pid), m_T(T), m_Q(Q) {
  if (!Utilities::isGoodAndPositive(Q)) throw std::runtime_error("secondary source vector is not valid.");
}

SecondarySource::~SecondarySource() { LOGD << "deleted SecondarySource for particle " << m_pid; }

double SecondarySource::get(const double& T) const {
  double value = 0;
  if (T > m_T.front() && T < m_T.back()) {
    value = GSL::LinearInterpolatorLog<double>(m_T, m_Q, T);
  }
  return value;
}

}  // namespace CRAMS