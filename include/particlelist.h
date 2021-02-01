#include <iomanip>
#include <map>
#include <vector>

#include "pid.h"

namespace CRAMS {

struct NucleusParameters {
  double abundance;
  double slope;
  double isotopicFractionISM;
  double decayTime;
  bool isStable;
  bool doPropagate;

  friend std::ostream& operator<<(std::ostream& stream, const NucleusParameters& inj) {
    stream << "(" << std::scientific << std::setprecision(3) << inj.abundance << ",";
    stream << std::fixed << std::setprecision(3) << inj.slope << ",";
    stream << std::boolalpha << inj.isStable << ")";
    return stream;
  }
};

typedef std::map<PID, NucleusParameters> List;

class ParticleList {
 private:
  List m_list;
  const std::string nucleilistFilename = "data/nucleilist.csv";

 public:
  ParticleList();
  virtual ~ParticleList();

  const List& getList() const { return m_list; }
  List& getList() { return m_list; }

  bool insert(const PID& key, const NucleusParameters& params);
  void setAbundance(const PID& key, const double& value);
  void setSlope(const PID& key, const double& value);
  void print() const;
  void readParamsFromFile(const std::string& filename);

 protected:
  void loadNucleilist(const std::string& filename);
  void setAbundanceChargeGroup(const int& charge, const double& abundance);
  void setSlopeChargeGroup(const int& charge, const double& slope);
  void setSlopeNuclei(const int& minCharge, const double& slope);

  // void set_from_file(const std::string& filename);
};

}  // namespace CRAMS