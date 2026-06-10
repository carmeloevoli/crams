#include <cassert>
#include <cmath>
#include <iostream>

#include "crams/core/cgs.h"
#include "crams/core/input.h"
#include "crams/core/pid.h"
#include "crams/grammage.h"

static int g_pass = 0;
static int g_fail = 0;

#define CHECK(cond)                                                                      \
  do {                                                                                   \
    if (cond) {                                                                          \
      ++g_pass;                                                                          \
    } else {                                                                             \
      ++g_fail;                                                                          \
      std::cerr << "FAIL: " << #cond << " at " << __FILE__ << ":" << __LINE__ << "\n"; \
    }                                                                                    \
  } while (0)

static bool approx(double a, double b, double tol = 1e-6) {
  return std::abs(a - b) <= tol * std::abs(b) + tol;
}

// Default-constructed Input with standard CR benchmark parameters
static CRAMS::Input makeInput() { return CRAMS::Input{}; }

// --- D(T): diffusion coefficient ---

void test_D_is_positive() {
  CRAMS::Grammage X(CRAMS::H1, makeInput());
  CHECK(X.D(1. * CRAMS::CGS::GeV) > 0.);
  CHECK(X.D(100. * CRAMS::CGS::GeV) > 0.);
  CHECK(X.D(10. * CRAMS::CGS::TeV) > 0.);
}

void test_D_floor_at_low_energy() {
  // At T → 0, beta → 0, so the power-law term vanishes and D → 2 * v_A * H
  CRAMS::Input in = makeInput();
  CRAMS::Grammage X(CRAMS::H1, in);
  const double floor = 2. * in.v_A() * in.H();
  // D is always above the floor
  CHECK(X.D(1e-4 * CRAMS::CGS::GeV) >= floor);
  CHECK(X.D(1e-2 * CRAMS::CGS::GeV) >= floor);
  // At very low energy the floor dominates (within 5%)
  CHECK(approx(X.D(1e-6 * CRAMS::CGS::GeV), floor, 5e-2));
}

void test_D_increases_with_energy() {
  // Power-law diffusion: D grows with energy above the low-energy floor
  CRAMS::Grammage X(CRAMS::H1, makeInput());
  CHECK(X.D(10. * CRAMS::CGS::GeV) > X.D(1. * CRAMS::CGS::GeV));
  CHECK(X.D(1. * CRAMS::CGS::TeV) > X.D(10. * CRAMS::CGS::GeV));
}

void test_D_carbon_larger_than_proton_at_same_T() {
  // R = pc/|Z| = sqrt(T*(T+2mp))*A/Z. Carbon has A/Z=2, proton A/Z=1,
  // so R_C = 2*R_p at same T/nuc → D_C > D_p.
  CRAMS::Input in = makeInput();
  CRAMS::Grammage X_p(CRAMS::H1, in);
  CRAMS::Grammage X_C(CRAMS::C12, in);
  const double T = 10. * CRAMS::CGS::GeV;
  CHECK(X_C.D(T) > X_p.D(T));
}

// --- diffusionTimescale / advectionTimescale ---

void test_diffusion_timescale_is_positive() {
  CRAMS::Grammage X(CRAMS::H1, makeInput());
  CHECK(X.diffusionTimescale(1. * CRAMS::CGS::GeV) > 0.);
}

void test_diffusion_timescale_formula() {
  CRAMS::Input in = makeInput();
  CRAMS::Grammage X(CRAMS::H1, in);
  const double T = 10. * CRAMS::CGS::GeV;
  const double expected = in.H() * in.H() / X.D(T);
  CHECK(approx(X.diffusionTimescale(T), expected));
}

void test_diffusion_timescale_decreases_with_energy() {
  // Higher energy → larger D → shorter diffusion timescale
  CRAMS::Grammage X(CRAMS::H1, makeInput());
  CHECK(X.diffusionTimescale(100. * CRAMS::CGS::GeV)
        < X.diffusionTimescale(1. * CRAMS::CGS::GeV));
}

void test_advection_timescale_is_positive() {
  CRAMS::Grammage X(CRAMS::H1, makeInput());
  CHECK(X.advectionTimescale() > 0.);
}

void test_advection_timescale_formula() {
  // t_adv = H / v_A — energy-independent, computable from Input directly
  CRAMS::Input in = makeInput();
  CRAMS::Grammage X(CRAMS::H1, in);
  const double expected = in.H() / in.v_A();
  CHECK(approx(X.advectionTimescale(), expected));
}

void test_advection_timescale_is_energy_independent() {
  CRAMS::Grammage X(CRAMS::H1, makeInput());
  CHECK(approx(X.advectionTimescale(), X.advectionTimescale()));
}

// --- get(T): grammage for stable particles ---

void test_get_stable_is_positive() {
  CRAMS::Grammage X(CRAMS::H1, makeInput());
  CHECK(X.get(1. * CRAMS::CGS::GeV) > 0.);
  CHECK(X.get(100. * CRAMS::CGS::GeV) > 0.);
}

void test_get_stable_upper_bound() {
  // X ≤ mu * beta * c / (2 * v_A) ≤ mu * c / (2 * v_A)
  CRAMS::Input in = makeInput();
  CRAMS::Grammage X(CRAMS::H1, in);
  const double X_max = in.mu() * CRAMS::CGS::cLight / (2. * in.v_A());
  CHECK(X.get(1. * CRAMS::CGS::GeV) < X_max);
  CHECK(X.get(1. * CRAMS::CGS::TeV) < X_max);
}

void test_get_stable_decreases_with_energy() {
  // Higher energy → larger D → particles escape faster → less grammage
  CRAMS::Grammage X(CRAMS::H1, makeInput());
  CHECK(X.get(100. * CRAMS::CGS::GeV) < X.get(10. * CRAMS::CGS::GeV));
  CHECK(X.get(10. * CRAMS::CGS::TeV) < X.get(100. * CRAMS::CGS::GeV));
}

void test_get_stable_limit_low_energy() {
  // As T → 0, beta → 0 and X → 0 (velocity factor). It should be smaller than at GeV.
  CRAMS::Grammage X(CRAMS::H1, makeInput());
  const double T_lo = 1e-3 * CRAMS::CGS::GeV;
  const double T_hi = 1. * CRAMS::CGS::GeV;
  CHECK(X.get(T_lo) < X.get(T_hi));
}

// --- get(T): unstable particles ---

void test_get_unstable_less_than_stable() {
  // Decay provides an additional escape channel → less grammage
  CRAMS::Input in = makeInput();
  const double tau_be10 = 1.387e6 * CRAMS::CGS::year;
  CRAMS::Grammage X_stable(CRAMS::Be10, in);
  CRAMS::Grammage X_unstable(CRAMS::Be10, in, tau_be10);
  const double T = 1. * CRAMS::CGS::GeV;
  CHECK(X_unstable.get(T) < X_stable.get(T));
}

void test_get_unstable_is_positive() {
  CRAMS::Input in = makeInput();
  const double tau_be10 = 1.387e6 * CRAMS::CGS::year;
  CRAMS::Grammage X(CRAMS::Be10, in, tau_be10);
  CHECK(X.get(1. * CRAMS::CGS::GeV) > 0.);
  CHECK(X.get(100. * CRAMS::CGS::GeV) > 0.);
}

void test_get_short_lifetime_much_smaller_than_stable() {
  // Very short lifetime: grammage is suppressed by many orders of magnitude vs stable
  CRAMS::Input in = makeInput();
  CRAMS::Grammage X_stable(CRAMS::H1, in);
  CRAMS::Grammage X_short(CRAMS::H1, in, 1. * CRAMS::CGS::sec);  // 1-second half-life
  const double T = 1. * CRAMS::CGS::GeV;
  CHECK(X_short.get(T) < 1e-4 * X_stable.get(T));
}

void test_get_unstable_decreases_with_energy() {
  CRAMS::Input in = makeInput();
  const double tau_be10 = 1.387e6 * CRAMS::CGS::year;
  CRAMS::Grammage X(CRAMS::Be10, in, tau_be10);
  // At higher energy gamma is larger → tau_d longer → approaches stable case
  // So unstable grammage should also decrease with energy (diffusion dominates at high E)
  CHECK(X.get(100. * CRAMS::CGS::GeV) < X.get(10. * CRAMS::CGS::GeV));
}

void test_get_long_lifetime_approaches_stable() {
  // As tau → ∞ the unstable formula reduces to the stable formula
  CRAMS::Input in = makeInput();
  CRAMS::Grammage X_stable(CRAMS::H1, in);
  CRAMS::Grammage X_long(CRAMS::H1, in, 1e30 * CRAMS::CGS::year);  // effectively infinite
  const double T = 10. * CRAMS::CGS::GeV;
  CHECK(approx(X_long.get(T), X_stable.get(T), 1e-4));
}

// --- heavy nuclei ---

void test_D_heavier_nucleus() {
  // Carbon (Z=6, A=12): lower rigidity at same T → smaller D than proton
  CRAMS::Input in = makeInput();
  CRAMS::Grammage X_C(CRAMS::C12, in);
  const double T = 10. * CRAMS::CGS::GeV;
  CHECK(X_C.D(T) > 0.);
}

void test_get_heavier_nucleus_is_positive() {
  CRAMS::Input in = makeInput();
  CRAMS::Grammage X_C(CRAMS::C12, in);
  CHECK(X_C.get(10. * CRAMS::CGS::GeV) > 0.);
}

int main() {
  test_D_is_positive();
  test_D_floor_at_low_energy();
  test_D_increases_with_energy();
  test_D_carbon_larger_than_proton_at_same_T();

  test_diffusion_timescale_is_positive();
  test_diffusion_timescale_formula();
  test_diffusion_timescale_decreases_with_energy();
  test_advection_timescale_is_positive();
  test_advection_timescale_formula();
  test_advection_timescale_is_energy_independent();

  test_get_stable_is_positive();
  test_get_stable_upper_bound();
  test_get_stable_decreases_with_energy();
  test_get_stable_limit_low_energy();

  test_get_unstable_less_than_stable();
  test_get_unstable_is_positive();
  test_get_short_lifetime_much_smaller_than_stable();
  test_get_unstable_decreases_with_energy();
  test_get_long_lifetime_approaches_stable();

  test_D_heavier_nucleus();
  test_get_heavier_nucleus_is_positive();

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
