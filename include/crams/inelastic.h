#ifndef CRAMS_INELASTIC_H_
#define CRAMS_INELASTIC_H_

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "crams/core/cgs.h"
#include "crams/core/pid.h"

namespace CRAMS {

double sigma_pp(const double& T);
double sigma_ST(const double& T, const int& A);

class InelasticXsec {
 public:
  virtual ~InelasticXsec();
  double getXsecOnISM(const PID& projectile, const double& T) const;
  virtual double getXsecOnHtarget(const PID& projectile, const double& T) const = 0;
};

class InXsecTripathi99 : public InelasticXsec {
 public:
  InXsecTripathi99();
  double getXsecOnHtarget(const PID& projectile, const double& T) const override;

 protected:
  void loadXsecTable(const std::string& filename);
  void buildEnergyArray();

 protected:
  std::map<PID, std::vector<double>> m_table;
  std::vector<double> m_T;
  const double m_T_min = 0.1 * CGS::GeV;
  const double m_T_max = 1e5 * CGS::GeV;
  const size_t m_T_size = 192;
  const std::string m_tableFilename = "data/crams_inelastic_tripathi1999.txt";
};

class InelasticXsecST98 : public InelasticXsec {
 public:
  InelasticXsecST98() = default;
  double getXsecOnHtarget(const PID& projectile, const double& T) const override;
};

}  // namespace CRAMS

#endif  // CRAMS_INELASTIC_H_
