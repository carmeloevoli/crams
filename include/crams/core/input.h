#ifndef CRAMS_CORE_INPUT_H_
#define CRAMS_CORE_INPUT_H_

#include <string>

#include "crams/core/cgs.h"
#include "crams/particlelist.h"

namespace CRAMS {

enum class FluxSolver {
  Analytical,
  CrankNicolson,
  Exponential,
};

class Input {
 public:
  Input() = default;
  ~Input() = default;

  void print() const;
  void readParamsFromFile(const std::string& filename);
  void setParam(const std::string& key, double value);
  void setSimname(const std::string& inifilename);

  double TSimMin() const { return m_TSimMin; }
  double TSimMax() const { return m_TSimMax; }
  size_t TSimSize() const { return m_TSimSize; }

  double ROutputMin() const { return m_ROutputMin; }
  double ROutputMax() const { return m_ROutputMax; }
  size_t ROutputSize() const { return m_ROutputSize; }

  bool doSecondary() const { return m_doSecondary; }

  double H() const { return m_H; }
  double mu() const { return m_mu; }
  double v_A() const { return m_v_A; }
  double R_b() const { return m_R_b; }
  double delta() const { return m_delta; }
  double ddelta() const { return m_ddelta; }
  double smoothness() const { return m_smoothness; }
  double D_0() const { return m_D_0; }
  double X_s() const { return m_X_s; }
  double modulationPotential() const { return m_modulationPotential; }
  double a_C() const { return m_a_C; }
  double a_D() const { return m_a_D; }
  size_t id() const { return m_id; }
  FluxSolver fluxSolver() const { return m_fluxSolver; }
  std::string fluxSolverName() const;
  const std::string& simname() const { return m_simname; }

 private:
  double m_TSimMin = 0.01 * CGS::GeV;
  double m_TSimMax = 10. * CGS::TeV;
  size_t m_TSimSize = 5 * 32 * 3;

  double m_ROutputMin = 0.1 * CGS::GeV;
  double m_ROutputMax = 10. * CGS::TeV;
  size_t m_ROutputSize = 100;

  bool m_doSecondary = true;

  double m_H = 7. * CGS::kpc;
  double m_mu = 2.3 * CGS::mgram / CGS::cm2;
  double m_v_A = 4.40940 * CGS::km / CGS::sec;
  double m_R_b = 290. * CGS::GeV;
  double m_delta = 5.65132e-01;
  double m_ddelta = 0.22;
  double m_smoothness = 0.1;
  double m_D_0 = 2.48255e28 * CGS::cm2 / CGS::sec;
  double m_X_s = -1.;
  double m_modulationPotential = 4.87754e-01 * CGS::GeV;
  double m_a_C = 1.0;
  double m_a_D = 1.0;
  size_t m_id = 0;
  FluxSolver m_fluxSolver = FluxSolver::CrankNicolson;
  std::string m_simname = "test";
};

}  // namespace CRAMS

#endif  // CRAMS_CORE_INPUT_H_
