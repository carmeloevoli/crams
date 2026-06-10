#include "crams/inelastic.h"

#include <plog/Log.h>

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <sstream>

#include "crams/utils/numeric.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

using Utilities::pow2;
using Utilities::pow3;

double sigma_pp(const double& T) {
  constexpr double E_threshold = 0.2797 * CGS::GeV;
  const double x = T / E_threshold;
  double value = 0;
  if (x > 1) {
    value = 30.7 - 0.96 * log(x) + 0.18 * pow2(log(x));
    value *= pow3(1 - pow(x, -1.9));
  }
  return value * CGS::mbarn;
}

double sigma_ST(const double& T, const int& A) {
  const double T_MeV = T / CGS::MeV;
  double value = 45. * std::pow((double)A, 0.7);
  value *= 1. + 0.016 * std::sin(5.3 - 2.63 * std::log(A));
  value *= 1. - 0.62 * std::exp(-T_MeV / 200.) * std::sin(10.9 * pow(T_MeV, -0.28));
  return value * CGS::mbarn;
}

InelasticXsec::~InelasticXsec() { LOGD << "deleted InelasticXsec"; }

double InelasticXsec::getXsecOnISM(const PID& projectile, const double& T) const {
  const double sigma_H = getXsecOnHtarget(projectile, T);
  return sigma_H * (1. + CGS::K_He * CGS::f_He) / (1. + CGS::f_He);
}

InXsecTripathi99::InXsecTripathi99() {
  buildEnergyArray();
  if (Utilities::fileExists(m_tableFilename))
    loadXsecTable(m_tableFilename);
  else
    throw std::runtime_error("Tripathi1999 inelastic xsecs file not found: " + m_tableFilename);
  if (m_table.empty()) throw std::runtime_error("Tripathi1999 inelastic xsecs table is empty: " + m_tableFilename);
  LOGD << "Tripathi1999 inelastic table read with " << m_table.size() << " projectiles";
}

void InXsecTripathi99::buildEnergyArray() {
  const double logRatio = std::log(m_T_max / m_T_min) / (m_T_size - 1);
  for (size_t i = 0; i < m_T_size; ++i) m_T.push_back(m_T_min * std::exp(i * logRatio));
}

double InXsecTripathi99::getXsecOnHtarget(const PID& projectile, const double& T) const {
  double sigma = 0;
  if (projectile.getZ() == 1 && projectile.getA() == 1) {
    sigma = sigma_pp(T);
  } else {
    const auto it = m_table.find(projectile);
    if (it == m_table.end())
      throw std::runtime_error("Tripathi1999 inelastic xsec not found for projectile " + projectile.toString());
    const double T_now = std::min(std::max(T, m_T.front()), m_T.back());
    sigma = Numeric::LinearInterpolator<double>(m_T, it->second, T_now);
  }
  return std::max(sigma, 1e-10 * CGS::mbarn);
}

void InXsecTripathi99::loadXsecTable(const std::string& filename) {
  std::ifstream inf(filename.c_str());
  std::string line;
  int Z_proj, A_proj;
  double x_temp;
  while (std::getline(inf, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    if (!(iss >> Z_proj >> A_proj)) continue;
    std::vector<double> x;
    x.reserve(m_T_size);
    for (size_t i = 0; i < m_T_size; ++i) {
      if (!(iss >> x_temp)) throw std::runtime_error("malformed Tripathi1999 inelastic xsecs row: " + line);
      x.emplace_back(x_temp * CGS::mbarn);
    }
    m_table[PID(Z_proj, A_proj)] = x;
  }
  inf.close();
}

}  // namespace CRAMS
