#ifndef CRAMS_RUNNER_H_
#define CRAMS_RUNNER_H_

#include "crams/core/cgs.h"
#include "crams/core/input.h"
#include "crams/core/output.h"
#include "crams/core/pid.h"
#include "crams/fragmentation.h"
#include "crams/inelastic.h"
#include "crams/particle.h"
#include "crams/particlelist.h"
#include "crams/utils/logging.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

class Runner {
 private:
  InelasticModel inelasticModel;
  FragmentationModel fragmentationModel;

  std::unique_ptr<InelasticXsec> inelasticXsecs;
  std::unique_ptr<NucFragXsec> nucfragXsecs;

  ParticleList injection;

 public:
  Runner(InelasticModel inelasticModel, FragmentationModel fragmentationModel, ParticleList injection);
  ~Runner() = default;
  void setInjectionParams(std::vector<double> abundances, std::vector<double> slopes);
  RigiditySpectra compute(Input input, bool dumpToFile = false, bool verbose = false, bool ignoreInputInitParams = false);
};

}  // namespace CRAMS

#endif  // CRAMS_RUNNER_H_
