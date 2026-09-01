#ifndef CRAMS_CORE_INPUT_H_
#define CRAMS_CORE_INPUT_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "crams/core/cgs.h"
#include "crams/core/pid.h"
#include "crams/particlelist.h"
#include "crams/physics/primary.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

struct SourceSpectrumFeature {
  virtual ~SourceSpectrumFeature() = default;
  virtual double rigidity() const = 0;
  virtual std::unique_ptr<SpectralFeature> toPrimarySourceFeature(const PID& pid) const = 0;
};

struct SourceSpectrumBreak : public SourceSpectrumFeature {
  SourceSpectrumBreak(double R_GV, double deltaSlope, double omega)
      : m_rigidity{R_GV * CGS::GeV}, m_deltaSlope{deltaSlope}, m_omega{omega} {}

  double rigidity() const override { return m_rigidity; }
  double deltaSlope() const { return m_deltaSlope; }
  double omega() const { return m_omega; }

  std::unique_ptr<SpectralFeature> toPrimarySourceFeature(const PID& pid) const override {
    const double T = Utilities::R2T(m_rigidity, pid);
    return std::make_unique<SpectralBreak>(T, m_deltaSlope, m_omega);
  }

 private:
  double m_rigidity;
  double m_deltaSlope;
  double m_omega;
};

struct SourceSpectrumLognormalFeature : public SourceSpectrumFeature {
  SourceSpectrumLognormalFeature(double R_GV, double sigma, double beta)
      : m_rigidity{R_GV * CGS::GeV}, m_sigma{sigma}, m_beta{beta} {}

  double rigidity() const override { return m_rigidity; }
  double sigma() const { return m_sigma; }
  double beta() const { return m_beta; }

  std::unique_ptr<SpectralFeature> toPrimarySourceFeature(const PID& pid) const override {
    const double T = Utilities::R2T(m_rigidity, pid);
    return std::make_unique<ErfcCutoff>(T, m_sigma, m_beta);
  }

 private:
  double m_rigidity;
  double m_sigma;
  double m_beta;
};

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
  void setTSim(double min_GeV, double max_GeV, int size) {
    m_TSimMin = min_GeV * CGS::GeV;
    m_TSimMax = max_GeV * CGS::GeV;
    m_TSimSize = size;
  }

  double ROutputMin() const { return m_ROutputMin; }
  double ROutputMax() const { return m_ROutputMax; }
  size_t ROutputSize() const { return m_ROutputSize; }
  void setROutput(double min_GV, double max_GV, int size) {
    m_ROutputMin = min_GV * CGS::GeV;
    m_ROutputMax = max_GV * CGS::GeV;
    m_ROutputSize = size;
  }

  void addSourceSpectrumFeature(std::shared_ptr<const SourceSpectrumFeature> feature) {
    m_sourceSpectrumFeatures.push_back(std::move(feature));
  }

  template <typename Feature, typename... Args>
  void addSourceSpectrumFeature(Args&&... args) {
    addSourceSpectrumFeature(std::make_shared<Feature>(std::forward<Args>(args)...));
  }

  const std::vector<std::shared_ptr<const SourceSpectrumFeature>>& sourceSpectrumFeatures() const {
    return m_sourceSpectrumFeatures;
  }

  // legacy setters for single-feature case

  void clearSourceSpectrumFeatures() { m_sourceSpectrumFeatures.clear(); }

  void setSourceSpectrumBreak(double R_GV, double deltaSlope, double omega) {
    clearSourceSpectrumFeatures();
    addSourceSpectrumFeature<SourceSpectrumBreak>(R_GV, deltaSlope, omega);
  }

  void setSourceSpectrumLognormal(double R_GV, double sigma, double beta) {
    clearSourceSpectrumFeatures();
    addSourceSpectrumFeature<SourceSpectrumLognormalFeature>(R_GV, sigma, beta);
  }

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

  // source spectral features, common for all primaries
  std::vector<std::shared_ptr<const SourceSpectrumFeature>> m_sourceSpectrumFeatures;

  bool m_doSecondary = true;

  // Multiplicative fudge on the Be isotope production cross-sections (1 = off).
  // Applied to the secondary source term in Particle::buildSecondarySource.
  double m_fudgeBe7 = 1.;
  double m_fudgeBe9 = 1.;
  double m_fudgeBe10 = 1.;

  FluxSolver m_fluxSolver = FluxSolver::CrankNicolson;
  InelasticModel m_inelasticModel = InelasticModel::Tripathi99;
  FragmentationModel m_fragmentationModel = FragmentationModel::UsineWebber03Coste12;

  size_t m_id = 0;
  std::string m_simname = "test";
};

}  // namespace CRAMS

#endif  // CRAMS_CORE_INPUT_H_
