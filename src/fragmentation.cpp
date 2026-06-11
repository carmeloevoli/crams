#include "crams/fragmentation.h"

#include <plog/Log.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "crams/utils/numeric.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

NucFragXsec::~NucFragXsec() { LOGD << "deleted NucFragXsec"; }

double NucFragXsec::getXsecOnISM(const PID& projectile, const PID& fragment, const double& T) const {
  const double sigma_H = getXsecOnHtarget(projectile, fragment, T);
  return sigma_H * (1. + CGS::K_He * CGS::f_He) / (1. + CGS::f_He);
}

NucFragFromTable::NucFragFromTable(std::string modelName, std::string tableFilename, double T_min, double T_max,
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
    throw std::runtime_error(m_modelName + " fragmentation xsecs file not found: " + m_tableFilename);
  if (m_table.empty())
    throw std::runtime_error(m_modelName + " fragmentation xsecs table is empty: " + m_tableFilename);
  LOGD << m_modelName << " fragmentation table read with " << m_table.size() << " channels";
}

void NucFragFromTable::buildEnergyArray() {
  const double logRatio = std::log(m_T_max / m_T_min) / (m_T_size - 1);
  for (size_t i = 0; i < m_T_size; ++i) m_T.push_back(m_T_min * std::exp(i * logRatio));
}

double NucFragFromTable::getXsecOnHtarget(const PID& projectile, const PID& fragment, const double& T) const {
  const auto it = m_table.find({projectile, fragment});
  if (it == m_table.end()) return 0.;  // channel not in this model -> no contribution
  // Small relative tolerance at the edges: the simulation energy grid can land
  // a hair outside [m_T_min, m_T_max] due to floating-point rounding in the
  // log-spaced grids, which is not a genuine out-of-range request.
  constexpr double edgeTol = 1e-6;
  if (T < m_T_min * (1. - edgeTol) || T > m_T_max * (1. + edgeTol))
    throw std::runtime_error(m_modelName + " fragmentation xsec requested at T = " + std::to_string(T / CGS::GeV) +
                             " GeV, outside the tabulated range [" + std::to_string(m_T_min / CGS::GeV) + ", " +
                             std::to_string(m_T_max / CGS::GeV) + "] GeV");
  // Clamp into the tabulated grid so a boundary point does not trip the interpolator.
  const double T_clamped = std::min(std::max(T, m_T.front()), m_T.back());
  return Numeric::LinearInterpolator<double>(m_T, it->second, T_clamped);
}

void NucFragFromTable::loadXsecTable() {
  std::ifstream inf(m_tableFilename.c_str());
  std::string line;
  int Z_proj, A_proj, Z_frag, A_frag;
  double x_temp;
  while (std::getline(inf, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    if (!(iss >> Z_proj >> A_proj >> Z_frag >> A_frag)) continue;
    std::vector<double> x;
    x.reserve(m_T_size);
    for (size_t i = 0; i < m_T_size; ++i) {
      if (!(iss >> x_temp)) throw std::runtime_error("malformed " + m_modelName + " fragmentation xsecs row: " + line);
      x.emplace_back(x_temp * CGS::mbarn);
    }
    m_table[{PID(Z_proj, A_proj), PID(Z_frag, A_frag)}] = x;
  }
  inf.close();
}

NucFragFluka4Dragon::NucFragFluka4Dragon()
    : NucFragFromTable("Fluka4Dragon", "data/crams_fragmentation_fluka4dragon.txt", 0.01 * CGS::GeV, 1e5 * CGS::GeV,
                       224) {}

}  // namespace CRAMS
