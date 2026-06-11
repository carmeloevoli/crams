#include "crams/inelastic.h"

#include <plog/Log.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "crams/utils/numeric.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

using Utilities::pow2;
using Utilities::pow3;

double sigma_pp(const double& T) {
  // Kafexhiu et al., Phys.Rev.D 90 (2014) 12, 123014
  // https://inspirehep.net/literature/1303850
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
  // R. Silberberg et al 1998 ApJ 501 911
  // https://iopscience.iop.org/article/10.1086/305862
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

InXsecFromTable::InXsecFromTable(std::string modelName, std::string tableFilename, double T_min, double T_max,
                                 size_t T_size)
    : m_modelName(std::move(modelName)),
      m_tableFilename(std::move(tableFilename)),
      m_T_min(T_min),
      m_T_max(T_max),
      m_T_size(T_size) {
  buildEnergyArray();
  if (Utilities::fileExists(m_tableFilename))
    loadXsecTable();
  else
    throw std::runtime_error(m_modelName + " inelastic xsecs file not found: " + m_tableFilename);
  if (m_table.empty()) throw std::runtime_error(m_modelName + " inelastic xsecs table is empty: " + m_tableFilename);
  LOGD << m_modelName << " inelastic table read with " << m_table.size() << " projectiles";
}

void InXsecFromTable::buildEnergyArray() {
  const double logRatio = std::log(m_T_max / m_T_min) / (m_T_size - 1);
  for (size_t i = 0; i < m_T_size; ++i) m_T.push_back(m_T_min * std::exp(i * logRatio));
}

double InXsecFromTable::getXsecOnHtarget(const PID& projectile, const double& T) const {
  double sigma = 0;
  if (projectile.getZ() == 1 && projectile.getA() == 1) {
    sigma = sigma_pp(T);
  } else {
    const auto it = m_table.find(projectile);
    if (it == m_table.end())
      throw std::runtime_error(m_modelName + " inelastic xsec not found for projectile " + projectile.toString());
    // Small relative tolerance at the edges: the simulation energy grid can land
    // a hair outside [m_T_min, m_T_max] due to floating-point rounding in the
    // log-spaced grids, which is not a genuine out-of-range request.
    constexpr double edgeTol = 1e-6;
    if (T < m_T_min * (1. - edgeTol) || T > m_T_max * (1. + edgeTol))
      throw std::runtime_error(m_modelName + " inelastic xsec requested at T = " + std::to_string(T / CGS::GeV) +
                               " GeV, outside the tabulated range [" + std::to_string(m_T_min / CGS::GeV) + ", " +
                               std::to_string(m_T_max / CGS::GeV) + "] GeV");
    // Clamp into the tabulated grid so a boundary point does not trip the interpolator.
    const double T_clamped = std::min(std::max(T, m_T.front()), m_T.back());
    sigma = Numeric::LinearInterpolator<double>(m_T, it->second, T_clamped);
  }
  return std::max(sigma, 1e-10 * CGS::mbarn);
}

void InXsecFromTable::loadXsecTable() {
  std::ifstream inf(m_tableFilename.c_str());
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
      if (!(iss >> x_temp)) throw std::runtime_error("malformed " + m_modelName + " inelastic xsecs row: " + line);
      x.emplace_back(x_temp * CGS::mbarn);
    }
    m_table[PID(Z_proj, A_proj)] = x;
  }
  inf.close();
}

InXsecTripathi99::InXsecTripathi99()
    : InXsecFromTable("Tripathi1999", "data/crams_inelastic_tripathi99.txt", 0.01 * CGS::GeV, 1e5 * CGS::GeV, 224) {}

InXsecGlauber::InXsecGlauber()
    : InXsecFromTable("Glauber", "data/crams_inelastic_glauber.txt", 0.01 * CGS::GeV, 1e5 * CGS::GeV, 224) {}

double InelasticXsecST98::getXsecOnHtarget(const PID& projectile, const double& T) const {
  if (projectile.getZ() == 1 && projectile.getA() == 1) return sigma_pp(T);
  return sigma_ST(T, projectile.getA());
}

}  // namespace CRAMS
