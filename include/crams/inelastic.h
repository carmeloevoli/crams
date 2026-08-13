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

// Inelastic cross-sections read from a tabulated XS4GCR file. Concrete models
// (Tripathi1999, Glauber, …) differ only in the filename and energy grid, so
// they share all the table machinery here and just supply those via the
// constructor.
class InXsecFromTable : public InelasticXsec {
 public:
  double getXsecOnHtarget(const PID& projectile, const double& T) const override;
  bool extrapolateToHighEnergies = true;

 protected:
  InXsecFromTable(std::string modelName, std::string tableFilename, double T_min, double T_max, size_t T_size);

 private:
  void buildEnergyArray();
  void loadXsecTable();

  const std::string m_modelName;
  const std::string m_tableFilename;
  const double m_T_min;
  const double m_T_max;
  const size_t m_T_size;
  std::map<PID, std::vector<double>> m_table;
  std::vector<double> m_T;
  std::vector<double> m_logT;
};

class InXsecTripathi99 : public InXsecFromTable {
 public:
  InXsecTripathi99();
};

class InXsecGlauber : public InXsecFromTable {
 public:
  InXsecGlauber();
};

class InXsecCrosec : public InXsecFromTable {
 public:
  InXsecCrosec();
};

class InelasticXsecST98 : public InelasticXsec {
 public:
  InelasticXsecST98() = default;
  double getXsecOnHtarget(const PID& projectile, const double& T) const override;
};

}  // namespace CRAMS

#endif  // CRAMS_INELASTIC_H_
