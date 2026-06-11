#ifndef CRAMS_PARTICLELIST_H_
#define CRAMS_PARTICLELIST_H_

#include <iomanip>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "crams/core/pid.h"

namespace CRAMS {

struct NucleusParameters {
  double abundance = 0.;
  double slope = 0.;
  double isotopicFractionISM = 0.;
  double decayTime = -1.;  // beta half-life at rest; negative sentinel means stable
  bool isStable = true;
  bool doPropagate = false;

  friend std::ostream& operator<<(std::ostream& stream, const NucleusParameters& inj) {
    stream << "(" << std::scientific << std::setprecision(3) << inj.abundance << ",";
    stream << std::fixed << std::setprecision(3) << inj.slope << ",";
    stream << std::boolalpha << inj.isStable << ")";
    return stream;
  }
};

using List = std::map<PID, NucleusParameters>;

class ParticleList {
 private:
  List m_list;
  std::map<int, std::vector<List::iterator>> m_particlesByCharge;
  bool m_rebuildChargeIndexBeforeUpdate = false;

 public:
  ParticleList();
  ParticleList(const ParticleList& other);
  ParticleList& operator=(const ParticleList& other);
  ParticleList(ParticleList&& other);
  ParticleList& operator=(ParticleList&& other);
  ~ParticleList() = default;

  const List& getList() const { return m_list; }
  List& getList() {
    m_rebuildChargeIndexBeforeUpdate = true;
    return m_list;
  }

  bool insert(const PID& key, const NucleusParameters& params);
  void setAbundance(const PID& key, double value);
  void setSlope(const PID& key, double value);
  void print() const;
  void readParamsFromFile(const std::string& filename);

 protected:
  void setParam(const std::string& key, double value);
  void loadNucleilist(const std::string& filename);
  void setAbundanceChargeGroup(int charge, double abundance);
  void setSlopeChargeGroup(int charge, double slope);
  void setSlopeNuclei(int minCharge, double slope);

 private:
  void applyDefaultInjectionParameters();
  void ensureChargeIndexFresh();
  void rebuildChargeIndex();
};

}  // namespace CRAMS

#endif  // CRAMS_PARTICLELIST_H_
