#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "crams/core/cgs.h"
#include "crams/core/input.h"
#include "crams/core/pid.h"
#include "crams/physics/primary.h"
#include "crams/utils/utilities.h"

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

#define CHECK_THROW(expr, exc) \
  do {                         \
    bool caught_ = false;      \
    try {                      \
      (void)(expr);            \
    } catch (const exc&) {     \
      caught_ = true;          \
    }                          \
    CHECK(caught_);            \
  } while (0)

static bool approx(double a, double b, double tol = 1e-6) { return std::abs(a - b) <= tol * std::abs(b) + tol; }

// Typical benchmark parameters
static const double slope = 4.5;
static const double abundance = 0.1;
static const double mu = CRAMS::Input{}.mu();

// --- zero abundance guard ---

void test_get_zero_abundance_returns_zero() {
  CRAMS::PrimarySource Q(CRAMS::H1, 0., slope, mu);
  CHECK(Q.get(1. * CRAMS::CGS::GeV) == 0.);
  CHECK(Q.get(1. * CRAMS::CGS::TeV) == 0.);
}

void test_get_negative_abundance_returns_zero() {
  CRAMS::PrimarySource Q(CRAMS::H1, -1., slope, mu);
  CHECK(Q.get(1. * CRAMS::CGS::GeV) == 0.);
}

// --- basic properties ---

void test_get_is_positive() {
  CRAMS::PrimarySource Q(CRAMS::H1, abundance, slope, mu);
  CHECK(Q.get(1. * CRAMS::CGS::GeV) > 0.);
  CHECK(Q.get(100. * CRAMS::CGS::GeV) > 0.);
  CHECK(Q.get(10. * CRAMS::CGS::TeV) > 0.);
}

void test_get_decreases_with_energy() {
  // Slope > 2, so 2-slope < 0: spectrum falls with momentum
  CRAMS::PrimarySource Q(CRAMS::H1, abundance, slope, mu);
  CHECK(Q.get(100. * CRAMS::CGS::GeV) < Q.get(10. * CRAMS::CGS::GeV));
  CHECK(Q.get(10. * CRAMS::CGS::TeV) < Q.get(100. * CRAMS::CGS::GeV));
}

// --- linearity in abundance ---

void test_get_scales_linearly_with_abundance() {
  const double T = 10. * CRAMS::CGS::GeV;
  CRAMS::PrimarySource Q1(CRAMS::H1, abundance, slope, mu);
  CRAMS::PrimarySource Q2(CRAMS::H1, 2. * abundance, slope, mu);
  CHECK(approx(Q2.get(T), 2. * Q1.get(T)));
}

void test_get_scales_linearly_with_abundance_at_high_energy() {
  const double T = 10. * CRAMS::CGS::TeV;
  CRAMS::PrimarySource Q1(CRAMS::H1, abundance, slope, mu);
  CRAMS::PrimarySource Q2(CRAMS::H1, 3. * abundance, slope, mu);
  CHECK(approx(Q2.get(T), 3. * Q1.get(T)));
}

// --- power-law spectrum ---

void test_get_power_law_at_ultrarelativistic_energy() {
  // At T >> mpc2: beta ≈ 1, pc ≈ T, so get(T1)/get(T2) ≈ (T1/T2)^(2-slope)
  // Using 10 TeV and 100 TeV: corrections are < 0.01%
  CRAMS::PrimarySource Q(CRAMS::H1, abundance, slope, mu);
  const double T1 = 10. * CRAMS::CGS::TeV;
  const double T2 = 100. * CRAMS::CGS::TeV;
  const double ratio = Q.get(T1) / Q.get(T2);
  // Expected: (T1/T2)^(2-slope) = (0.1)^(-2.5) = 10^2.5 = 316.228
  const double expected = std::pow(T1 / T2, 2. - slope);
  CHECK(approx(ratio, expected, 1e-3));
}

void test_get_power_law_ratio_matches_momentum_ratio() {
  // More general: get(T1)/get(T2) = (beta2/beta1) * (pc1/pc2)^(2-slope)
  CRAMS::PrimarySource Q(CRAMS::H1, abundance, slope, mu);
  const double T1 = 1. * CRAMS::CGS::GeV;
  const double T2 = 10. * CRAMS::CGS::GeV;
  const double ratio = Q.get(T1) / Q.get(T2);
  const double pc1 = CRAMS::Utilities::T2pc(T1, CRAMS::H1);
  const double pc2 = CRAMS::Utilities::T2pc(T2, CRAMS::H1);
  const double beta1 = CRAMS::Utilities::T2beta(T1);
  const double beta2 = CRAMS::Utilities::T2beta(T2);
  const double expected = (beta2 / beta1) * std::pow(pc1 / pc2, 2. - slope);
  CHECK(approx(ratio, expected, 1e-9));
}

void test_get_steeper_slope_falls_faster() {
  // Steeper slope → faster energy fall-off
  const double T1 = 1. * CRAMS::CGS::GeV;
  const double T2 = 100. * CRAMS::CGS::GeV;
  CRAMS::PrimarySource Q_soft(CRAMS::H1, abundance, 4.2, mu);
  CRAMS::PrimarySource Q_hard(CRAMS::H1, abundance, 4.8, mu);
  // ratio at (T2/T1) for steeper slope must be smaller (more suppressed)
  const double r_soft = Q_soft.get(T2) / Q_soft.get(T1);
  const double r_hard = Q_hard.get(T2) / Q_hard.get(T1);
  CHECK(r_hard < r_soft);
}

void test_get_combines_multiple_features() {
  CRAMS::PrimarySource Q(CRAMS::H1, abundance, slope, mu);
  const double T_break = 10. * CRAMS::CGS::GeV;
  const double T_cut = 100. * CRAMS::CGS::GeV;

  Q.addFeature(std::make_unique<CRAMS::SpectralBreak>(T_break, 1.5, 0.5));
  Q.addFeature(std::make_unique<CRAMS::ErfcCutoff>(T_cut, 0.3, 0.7));

  const double base = Q.get(1. * CRAMS::CGS::GeV);
  const double with_break = Q.get(10. * T_break);
  const double with_both = Q.get(10. * T_cut);

  CHECK(base > 0.);
  CHECK(with_break < base);
  CHECK(with_both < with_break);
}

// --- mass-number scaling ---

void test_get_A_scaling() {
  // m_norm ∝ A, pc ∝ A → get ∝ A * A^(2-slope) = A^(3-slope)
  // For slope=4.5: get ∝ A^(-1.5), so proton > carbon at same T/nuc and equal abundance
  const double T = 10. * CRAMS::CGS::GeV;
  CRAMS::PrimarySource Q_p(CRAMS::H1, abundance, slope, mu);
  CRAMS::PrimarySource Q_C(CRAMS::C12, abundance, slope, mu);
  // Carbon A=12: get_C ∝ 12^(3-4.5) = 12^(-1.5) = 1/41.6 smaller factor
  const double expected_ratio = std::pow(12., 3. - slope);  // < 1 for slope > 3
  CHECK(approx(Q_C.get(T) / Q_p.get(T), expected_ratio, 1e-6));
}

void test_get_heavier_nucleus_lower_for_steep_slope() {
  // For typical slope > 3, heavier nuclei have lower get at same T/nuc and equal abundance
  const double T = 10. * CRAMS::CGS::GeV;
  CRAMS::PrimarySource Q_p(CRAMS::H1, abundance, slope, mu);
  CRAMS::PrimarySource Q_C(CRAMS::C12, abundance, slope, mu);
  CHECK(Q_p.get(T) > Q_C.get(T));
}

// --- gas surface density scaling ---

void test_get_scales_inversely_with_surface_density() {
  // m_norm ∝ 1/surfaceDensity → doubling mu halves get(T)
  const double T = 10. * CRAMS::CGS::GeV;
  CRAMS::PrimarySource Q1(CRAMS::H1, abundance, slope, mu);
  CRAMS::PrimarySource Q2(CRAMS::H1, abundance, slope, 2. * mu);
  CHECK(approx(Q2.get(T), 0.5 * Q1.get(T)));
}

// --- slope validity ---

void test_invalid_slope_too_low_throws() {
  // GammaIntegral requires slope in (4, 5); abundance > 0 triggers the call
  CHECK_THROW(CRAMS::PrimarySource(CRAMS::H1, abundance, 3.5, mu), std::invalid_argument);
}

void test_invalid_slope_too_high_throws() {
  CHECK_THROW(CRAMS::PrimarySource(CRAMS::H1, abundance, 5.5, mu), std::invalid_argument);
}

void test_invalid_slope_boundary_low_throws() {
  // slope == 4.0 is not in the open interval (4, 5)
  CHECK_THROW(CRAMS::PrimarySource(CRAMS::H1, abundance, 4.0, mu), std::invalid_argument);
}

void test_invalid_slope_ignored_when_zero_abundance() {
  // Zero abundance skips GammaIntegral, so no throw even for invalid slope
  CRAMS::PrimarySource Q(CRAMS::H1, 0., 3.0, mu);
  CHECK(Q.get(1. * CRAMS::CGS::GeV) == 0.);
}

int main() {
  test_get_zero_abundance_returns_zero();
  test_get_negative_abundance_returns_zero();

  test_get_is_positive();
  test_get_decreases_with_energy();

  test_get_scales_linearly_with_abundance();
  test_get_scales_linearly_with_abundance_at_high_energy();

  test_get_power_law_at_ultrarelativistic_energy();
  test_get_power_law_ratio_matches_momentum_ratio();
  test_get_steeper_slope_falls_faster();
  test_get_combines_multiple_features();

  test_get_A_scaling();
  test_get_heavier_nucleus_lower_for_steep_slope();
  test_get_scales_inversely_with_surface_density();

  test_invalid_slope_too_low_throws();
  test_invalid_slope_too_high_throws();
  test_invalid_slope_boundary_low_throws();
  test_invalid_slope_ignored_when_zero_abundance();

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
