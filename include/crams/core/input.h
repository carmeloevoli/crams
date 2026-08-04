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

enum class InelasticModel {
  Tripathi99,
  Glauber,
  Crosec,
};

enum class FragmentationModel {
  Fluka4Dragon,
  UsineGalprop17Opt12,
  UsineGalprop17Opt22,
  UsineWebber03Coste12,
  Evoli2019,
  Evoli2026W93,
  Evoli2026St99,
};

InelasticModel parseInelasticModel(const std::string& value);
FragmentationModel parseFragmentationModel(const std::string& value);
FluxSolver parseFluxSolver(const std::string& value);

class Input {
 public:
  Input(double H_kpc = 7., double mu_mg_per_cm2 = 2.3, double v_A_km_sec = 4.40940, double R_b_GV = 290.,
        double delta = 5.65132e-01, double ddelta = 0.22, double smoothness = 0.1, double D_0_cm2_sec = 2.48255e28,
        double X_s = -1., double modulationPotential = 4.87754e-01)
      : m_H{H_kpc * CGS::kpc},
        m_mu{mu_mg_per_cm2 * CGS::mgram / CGS::cm2},
        m_v_A{v_A_km_sec * CGS::km / CGS::sec},
        m_R_b{R_b_GV * CGS::GeV},
        m_delta{delta},
        m_ddelta{ddelta},
        m_smoothness{smoothness},
        m_D_0{D_0_cm2_sec * CGS::cm2 / CGS::sec},
        m_X_s{X_s},
        m_modulationPotential{modulationPotential * CGS::GeV} {}

  ~Input() = default;

  std::string describe() const;
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
  double fudgeBe7() const { return m_fudgeBe7; }
  double fudgeBe9() const { return m_fudgeBe9; }
  double fudgeBe10() const { return m_fudgeBe10; }
  size_t id() const { return m_id; }
  FluxSolver fluxSolver() const { return m_fluxSolver; }
  std::string fluxSolverName() const;
  InelasticModel inelasticModel() const { return m_inelasticModel; }
  std::string inelasticModelName() const;
  FragmentationModel fragmentationModel() const { return m_fragmentationModel; }
  std::string fragmentationModelName() const;
  const std::string& simname() const { return m_simname; }

 private:
  // main physics params
  double m_H;
  double m_mu;
  double m_v_A;
  double m_R_b;
  double m_delta;
  double m_ddelta;
  double m_smoothness;
  double m_D_0;
  double m_X_s;
  double m_modulationPotential;

  // computation grid
  double m_TSimMin = 0.1 * CGS::GeV;
  double m_TSimMax = 100. * CGS::TeV;
  size_t m_TSimSize = 300;

  // output grid
  double m_ROutputMin = 1. * CGS::GeV;
  double m_ROutputMax = 10. * CGS::TeV;
  size_t m_ROutputSize = 100;

  bool m_doSecondary = true;

  // Multiplicative fudge on the Be isotope production cross-sections (1 = off).
  // Applied to the secondary source term in Particle::buildSecondarySource.
  double m_fudgeBe7 = 1.;
  double m_fudgeBe9 = 1.;
  double m_fudgeBe10 = 1.;
  size_t m_id = 0;
  FluxSolver m_fluxSolver = FluxSolver::CrankNicolson;
  InelasticModel m_inelasticModel = InelasticModel::Tripathi99;
  FragmentationModel m_fragmentationModel = FragmentationModel::UsineWebber03Coste12;
  std::string m_simname = "test";
};

}  // namespace CRAMS

#endif  // CRAMS_CORE_INPUT_H_
