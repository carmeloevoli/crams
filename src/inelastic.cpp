#include "inelastic.h"

#include <plog/Log.h>

#include <random>

#include "utilities.h"
#include "xsecs/Tripathi99.h"

namespace CRAMS {

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

InelasticXsec::InelasticXsec() {}

InelasticXsec::InelasticXsec(const PID& proj, bool doRandom) : m_proj(proj), m_doRandom(doRandom) {
  if (proj != H1) {
    m_randomFactor = (m_doRandom) ? Utilities::computeRandomFactor(m_randomFactorVariance) : 1.;
  }
}

InelasticXsec::~InelasticXsec() { LOGD << "deleted InelasticXsec for particle " << m_proj; }

double InelasticXsec::getXsecOnISM(const double& T) const {
  const double sigma_H = getXsecOnHtarget(T);
  return m_randomFactor * sigma_H * (1. + CGS::K_He * CGS::f_He) / (1. + CGS::f_He);
}

InXsecTripathi99::InXsecTripathi99() {}

InXsecTripathi99::InXsecTripathi99(const PID& proj, bool doError) : InelasticXsec(proj, doError) {}

double InXsecTripathi99::getXsecOnHtarget(const double& T) const {
  double sigma = 0;
  if (m_proj == H1)
    sigma = sigma_pp(T);
  else
    sigma = Tripathi99::inelastic_sigma(1, 1, m_proj.getA(), m_proj.getZ(), T);
  return std::max(sigma, 1e-10 * CGS::mbarn);
}

InXsecCROSEC::InXsecCROSEC() {}

InXsecCROSEC::InXsecCROSEC(const PID& proj, bool doError) : InelasticXsec(proj, doError) {
  buildEnergyArray();
  if (Utilities::fileExists(m_tableFilename))
    loadXsecTable(m_tableFilename);
  else
    throw std::runtime_error("CROSEC xsecs file not found");
}

void InXsecCROSEC::buildEnergyArray() {
  double T = m_T_min;
  for (size_t i = 0; i < m_T_size; ++i) {
    m_T.push_back(T);
    T *= m_T_ratio;
  }
}

double InXsecCROSEC::getXsecOnHtarget(const double& T) const {
  double sigma = 0;
  if (m_proj == H1)
    sigma = sigma_pp(T);
  else
    sigma = (T >= m_T.back()) ? m_table.back() : Utilities::LinearInterpolator(m_T, m_table, T);
  return std::max(sigma, 1e-10 * CGS::mbarn);
}

void InXsecCROSEC::loadXsecTable(const std::string& filename) {
  std::ifstream inf(filename.c_str());
  int Z_proj, A_proj;
  double x_temp;
  while (inf) {
    inf >> Z_proj >> A_proj;
    std::vector<double> x;
    x.reserve(m_T_size);
    for (size_t i = 0; i < m_T_size; ++i) {
      inf >> x_temp;
      x.emplace_back(x_temp * CGS::mbarn);
    }
    if (PID(Z_proj, A_proj) == m_proj && inf.good()) copy(x.begin(), x.end(), back_inserter(m_table));
  }
  inf.close();
  LOGD << "inelastic table read with " << m_table.size() << " points";
}

}  // namespace CRAMS