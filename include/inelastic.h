#ifndef INCLUDE_INELASTIC_H_
#define INCLUDE_INELASTIC_H_

#include <cmath>
#include <vector>

#include "cgs.h"
#include "pid.h"

namespace CRAMS {

double sigma_pp(const double& T);
double sigma_ST(const double& T);

class InelasticXsec {
 public:
  InelasticXsec();
  InelasticXsec(const PID& proj, bool doRandom = false);
  virtual ~InelasticXsec();
  double getXsecOnISM(const double& T) const;
  virtual double getXsecOnHtarget(const double& T) const = 0;

 protected:
  PID m_proj;
  double m_randomFactor = 1.;
  double m_randomFactorVariance = 0.05;
  bool m_doRandom = false;

 private:
  double getRandomFactor() const;
};

class InXsecTripathi99 : public InelasticXsec {
 public:
  InXsecTripathi99();
  InXsecTripathi99(const PID& proj, bool doRandom = false);
  double getXsecOnHtarget(const double& T) const override;
};

class InXsecCROSEC : public InelasticXsec {
 public:
  InXsecCROSEC();
  InXsecCROSEC(const PID& proj, bool doRandom = false);
  double getXsecOnHtarget(const double& T) const override;

 protected:
  void loadXsecTable(const std::string& filename);
  void buildEnergyArray();

 protected:
  std::vector<double> m_table;
  std::vector<double> m_T;
  const double m_T_min = 0.1 * CGS::GeV;
  const size_t m_T_size = 100;
  const double m_T_ratio = 1.1;
  const std::string m_tableFilename = "data/xsecs_total_CROSEC_0.1_100_1.1.txt";
};

}  // namespace CRAMS

#endif  // INCLUDE_INELASTIC_H_
