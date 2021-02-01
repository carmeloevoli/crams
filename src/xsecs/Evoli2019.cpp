#include "xsecs/Evoli2019.h"

#include <plog/Log.h>

#include "csvreader.h"
#include "utilities.h"

namespace CRAMS {

SpallationXsecs::SpallationXsecs(const PID& fragment, const double& fudgeFactor, bool doRandom)
    : m_fragment(fragment), m_fudgeFactor(fudgeFactor), m_doRandom(doRandom) {
  buildEnergyArray();
  if (Utilities::fileExists(m_tableFilename))
    loadXsecTable(m_tableFilename);
  else
    throw std::runtime_error("spallation xsecs file not found");
  LOGW << "read " << m_table.size() << " cross-sections for " << m_fragment;
}

SpallationXsecs::~SpallationXsecs() { LOGD << "deleted SpallationXsecs for fragment " << m_fragment; }

void SpallationXsecs::buildEnergyArray() {
  double T = m_T_min;
  for (size_t i = 0; i < m_T_size; ++i) {
    m_T.push_back(T);
    T *= m_T_ratio;
  }
}

double SpallationXsecs::getXsecOnISM(const PID& projectile, const double& T) const {
  double sigma_H = getXsecOnHtarget(projectile, T);
  if (m_fragment.getZ() == 4) sigma_H *= m_fudgeFactor;
  return sigma_H * (1. + CGS::K_He * CGS::f_He) / (1. + CGS::f_He);
}

double SpallationXsecs::getXsecOnHtarget(const PID& projectile, const double& T) const {
  double value = 0;
  auto it = m_table.find(projectile);
  if (it != m_table.end()) {
    double T_now = std::min(T, m_T.back());
    value = Utilities::LinearInterpolatorLog(m_T, it->second, T_now);
    auto it_error = m_randomFactors.find(projectile);
    value *= (m_doRandom) ? it_error->second : 1.;
  }
  return value;
}

void SpallationXsecs::loadXsecTable(const std::string& filename) {
  CSVReader reader(filename);
  auto xsecslist = reader.getData();
  for (const auto& xsec : xsecslist) {
    const int Z_frag = atoi(xsec[0].c_str());
    const int A_frag = atoi(xsec[1].c_str());
    const int Z_proj = atoi(xsec[2].c_str());
    const int A_proj = atoi(xsec[3].c_str());
    if (PID(Z_frag, A_frag) == m_fragment) {
      std::vector<double> x;
      x.reserve(m_T_size);
      for (size_t i = 0; i < m_T_size; ++i) {
        const auto x_temp = atof(xsec[4 + i].c_str());
        x.emplace_back(std::max(x_temp, 1e-10) * CGS::mbarn);
      }
      m_table[PID(Z_proj, A_proj)] = x;
      m_randomFactors[PID(Z_proj, A_proj)] = Utilities::computeRandomFactor(m_randomFactorVariance);
    }
  }
}

}  // namespace CRAMS
