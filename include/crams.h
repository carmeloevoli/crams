#ifndef CRAMS_H_
#define CRAMS_H_

// Physical constants and units
#include "crams/core/cgs.h"

// Core types
#include "crams/core/pid.h"
#include "crams/particlelist.h"

// Simulation components
#include "crams/core/input.h"
#include "crams/core/output.h"
#include "crams/fragmentation.h"
#include "crams/inelastic.h"
#include "crams/particle.h"

// Infrastructure
#include "crams/utils/logging.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

class Runner {
 private:
  std::unique_ptr<InXsecGlauber> inelasticXsecs;
  std::unique_ptr<NucFragXsec> nucfragXsecs;
  FluxSolver solver;

 public:
  Runner(FluxSolver solver, InelasticModel inelasticModel, FragmentationModel fragmentationModel);
  ~Runner() = default;
//   compute(init)
};

}  // namespace CRAMS

#endif  // CRAMS_H_
