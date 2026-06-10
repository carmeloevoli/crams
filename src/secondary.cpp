#include "crams/secondary.h"

#include <plog/Log.h>

#include <stdexcept>

#include "crams/utils/numeric.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

SecondarySource::SecondarySource(const PID& pid, const std::vector<double>& T,
                                 const std::vector<double>& Q)
    : m_pid(pid), m_T(T), m_Q(Q) {
  if (!Utilities::isGoodAndPositive(Q))
    throw std::runtime_error("secondary source vector is not valid");
}

SecondarySource::~SecondarySource() { LOGD << "deleted SecondarySource for particle " << m_pid; }

double SecondarySource::get(double T) const {
  if (T <= m_T.front() || T >= m_T.back()) return 0.;
  return Numeric::LinearInterpolatorLog<double>(m_T, m_Q, T);
}

}  // namespace CRAMS
