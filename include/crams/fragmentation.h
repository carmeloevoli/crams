#ifndef CRAMS_FRAGMENTATION_H_
#define CRAMS_FRAGMENTATION_H_

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "crams/core/cgs.h"
#include "crams/core/pid.h"

namespace CRAMS {

// Nuclear fragmentation (spallation) cross-sections for the channel
// projectile -> fragment. Mirrors the InelasticXsec hierarchy.
class NucFragXsec {
 public:
  virtual ~NucFragXsec();
  // Cross-section for projectile -> fragment on the ISM (H + He), at kinetic
  // energy per nucleon T. Returns 0 if the channel is absent from the model.
  double getXsecOnISM(const PID& projectile, const PID& fragment, const double& T) const;
  virtual double getXsecOnHtarget(const PID& projectile, const PID& fragment, const double& T) const = 0;
};

// Fragmentation cross-sections read from a tabulated XS4GCR file. Concrete
// models (Fluka4Dragon, …) differ only in the filename and energy grid, so they
// share all the table machinery here and just supply those via the constructor.
class NucFragFromTable : public NucFragXsec {
 public:
  double getXsecOnHtarget(const PID& projectile, const PID& fragment, const double& T) const override;

 protected:
  NucFragFromTable(std::string modelName, std::string tableFilename, double T_min, double T_max, size_t T_size);

 private:
  void buildEnergyArray();
  void loadXsecTable();

  const std::string m_modelName;
  const std::string m_tableFilename;
  const double m_T_min;
  const double m_T_max;
  const size_t m_T_size;
  std::map<std::pair<PID, PID>, std::vector<double>> m_table;  // (projectile, fragment) -> sigma(T)
  std::vector<double> m_T;
};

class NucFragFluka4Dragon : public NucFragFromTable {
 public:
  NucFragFluka4Dragon();
};

}  // namespace CRAMS

#endif  // CRAMS_FRAGMENTATION_H_
