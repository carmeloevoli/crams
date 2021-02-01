#ifndef INCLUDE_SPALLATION_H_
#define INCLUDE_SPALLATION_H_

#include <map>
#include <string>
#include <vector>

#include "cgs.h"
#include "pid.h"

namespace CRAMS {

class SpallationXsecs {
 public:
  SpallationXsecs(const PID& fragment, const double& fudgeFactor, bool doRandom = false);
  virtual ~SpallationXsecs();
  double getXsecOnISM(const PID& projectile, const double& T) const;

 protected:
  void buildEnergyArray();
  void loadXsecTable(const std::string& filename);
  double getXsecOnHtarget(const PID& projectile, const double& T) const;

 protected:
  bool m_doRandom = false;
  double m_randomFactor = 1;
  double m_randomFactorVariance = 0.3;
  PID m_fragment;
  std::map<PID, double> m_randomFactors;
  std::map<PID, std::vector<double> > m_table;
  std::vector<double> m_T;
  double m_fudgeFactor = 1;

#ifdef EVOLI2019
  std::string m_tableFilename = "data/crxsecs_fragmentation_Evoli2019_cumulative.csv";
  const double m_T_min = 0.01 * CGS::GeV;
  const size_t m_T_size = 160;
  const double m_T_ratio = 1.0751;
#endif
};

}  // namespace CRAMS

#endif  // INCLUDE_SPALLATION_H_