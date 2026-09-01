// Parses examples/crams.ini and checks that every documented key is understood
// and produces the expected value. This makes the example file a live
// specification: if a key is renamed/removed, or a value in the file is edited,
// this test fails. The .ini path is injected by CMake as EXAMPLE_INI.
#include <cmath>
#include <iostream>
#include <string>

#include "crams/core/cgs.h"
#include "crams/core/input.h"
#include "crams/particlelist.h"

#ifndef EXAMPLE_INI
#error "EXAMPLE_INI must be defined by the build system (path to examples/crams.ini)"
#endif

static int g_pass = 0;
static int g_fail = 0;

#define CHECK(cond)                                                                    \
  do {                                                                                 \
    if (cond) {                                                                        \
      ++g_pass;                                                                        \
    } else {                                                                           \
      ++g_fail;                                                                        \
      std::cerr << "FAIL: " << #cond << " at " << __FILE__ << ":" << __LINE__ << "\n"; \
    }                                                                                  \
  } while (0)

static bool approx(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol * std::abs(b) + tol; }

// Slope of the (Z, A) primary isotope, or -1 if not found.
static double slopeOf(const CRAMS::ParticleList& pl, int Z, int A) {
  for (const auto& entry : pl.getList())
    if (entry.first.getZ() == Z && entry.first.getA() == A && !entry.first.isTertiary()) return entry.second.slope;
  return -1.;
}

// Total source abundance summed over the isotopes of charge Z.
static double abundanceForZ(const CRAMS::ParticleList& pl, int Z) {
  double sum = 0.;
  for (const auto& entry : pl.getList())
    if (entry.first.getZ() == Z) sum += entry.second.abundance;
  return sum;
}

void test_input_transport() {
  CRAMS::Input in;
  in.readParamsFromFile(EXAMPLE_INI);
  CHECK(approx(in.D_0(), 3.0 * 1e28 * CRAMS::CGS::cm2 / CRAMS::CGS::sec));
  CHECK(approx(in.delta(), 0.50));
  CHECK(approx(in.ddelta(), 0.20));
  CHECK(approx(in.R_b(), 300.0 * CRAMS::CGS::GeV));
  CHECK(approx(in.v_A(), 5.0 * CRAMS::CGS::km / CRAMS::CGS::sec));
  CHECK(approx(in.H(), 5.0 * CRAMS::CGS::kpc));
  CHECK(approx(in.modulationPotential(), 0.60 * CRAMS::CGS::GeV));
  CHECK(in.X_s() <= 0.);  // xs = -1 disables source grammage
}

void test_input_numerics() {
  CRAMS::Input in;
  in.readParamsFromFile(EXAMPLE_INI);
  CHECK(in.id() == 7);
  CHECK(in.fluxSolver() == CRAMS::FluxSolver::CrankNicolson);
  CHECK(in.fluxSolverName() == "crank_nicolson");
  CHECK(in.inelasticModel() == CRAMS::InelasticModel::Glauber);
  CHECK(in.inelasticModelName() == "glauber");
  CHECK(in.fragmentationModel() == CRAMS::FragmentationModel::Fluka4Dragon);
  CHECK(in.fragmentationModelName() == "fluka4dragon");
}

void test_injection_slopes() {
  CRAMS::ParticleList pl;
  pl.readParamsFromFile(EXAMPLE_INI);
  CHECK(approx(slopeOf(pl, 1, 1), 4.40));   // hslope -> protons
  CHECK(approx(slopeOf(pl, 2, 4), 4.35));   // heslope -> helium
  CHECK(approx(slopeOf(pl, 6, 12), 4.30));  // slope   -> nuclei (Z >= 3), e.g. C12
  CHECK(approx(slopeOf(pl, 8, 16), 4.30));  // slope   -> nuclei (Z >= 3), e.g. O16
}

void test_injection_abundances() {
  CRAMS::ParticleList pl;
  pl.readParamsFromFile(EXAMPLE_INI);
  CHECK(abundanceForZ(pl, 1) > 0.);   // qh
  CHECK(abundanceForZ(pl, 6) > 0.);   // qc
  CHECK(abundanceForZ(pl, 26) > 0.);  // qfe
  CHECK(abundanceForZ(pl, 3) > 0.);   // qli set non-zero in the example (default 0)
}

void test_injection_break() {
  CRAMS::Input in;
  in.readParamsFromFile(EXAMPLE_INI);
  const auto& features = in.sourceSpectrumFeatures();
  CHECK(features.size() == 1);
  const auto* sourceBreak = dynamic_cast<const CRAMS::SourceSpectrumBreak*>(features[0].get());
  CHECK(sourceBreak != nullptr);
  CHECK(approx(sourceBreak->rigidity(), 13. * CRAMS::CGS::TeV));
  CHECK(approx(sourceBreak->deltaSlope(), 0.32));
  CHECK(approx(sourceBreak->omega(), 0.1));
}

int main() {
  test_input_transport();
  test_input_numerics();
  test_injection_slopes();
  test_injection_abundances();

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
