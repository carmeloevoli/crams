#include <cmath>
#include <iostream>

#include "crams/core/cgs.h"
#include "crams/core/input.h"
#include "crams/core/pid.h"
#include "crams/physics/losses.h"
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

static bool approx(double a, double b, double tol = 1e-6) { return std::abs(a - b) <= tol * std::abs(b) + tol; }

static CRAMS::Input makeInput() { return CRAMS::Input{}; }

// --- sign: all loss terms must be negative ---

void test_get_is_negative() {
  CRAMS::Losses L(CRAMS::H1, makeInput());
  CHECK(L.get(1. * CRAMS::CGS::GeV) < 0.);
  CHECK(L.get(10. * CRAMS::CGS::TeV) < 0.);
}

void test_dEdX_adiabatic_is_negative() {
  CRAMS::Losses L(CRAMS::H1, makeInput());
  CHECK(L.dEdX_adiabatic(1. * CRAMS::CGS::GeV) < 0.);
  CHECK(L.dEdX_adiabatic(1. * CRAMS::CGS::TeV) < 0.);
}

void test_dEdX_ionization_is_negative() {
  CRAMS::Losses L(CRAMS::H1, makeInput());
  CHECK(L.dEdX_ionization(1. * CRAMS::CGS::GeV) < 0.);
  CHECK(L.dEdX_ionization(1. * CRAMS::CGS::TeV) < 0.);
}

void test_dTdt_ionization_is_positive() {
  CRAMS::Losses L(CRAMS::H1, makeInput());
  CHECK(L.dTdt_ionization(1. * CRAMS::CGS::GeV, 1. / CRAMS::CGS::cm3) > 0.);
}

void test_get_equals_sum_of_components() {
  CRAMS::Losses L(CRAMS::H1, makeInput());
  const double T = 5. * CRAMS::CGS::GeV;
  CHECK(approx(L.get(T), L.dEdX_adiabatic(T) + L.dEdX_ionization(T)));
}

// --- dEdX_adiabatic: exact ratio test ---

void test_adiabatic_ratio_equals_momentum_ratio() {
  // |dEdX_adiabatic(T)| = factorAdv * pc, so the ratio must equal pc(T1)/pc(T2)
  CRAMS::Losses L(CRAMS::H1, makeInput());
  const double T1 = 1. * CRAMS::CGS::GeV;
  const double T2 = 100. * CRAMS::CGS::GeV;
  const double ratio = L.dEdX_adiabatic(T1) / L.dEdX_adiabatic(T2);
  // pc = sqrt(T*(T+2*mpc2)), same as T2pc with A=1
  const double pc1 = CRAMS::Utilities::T2pc(T1, CRAMS::H1);
  const double pc2 = CRAMS::Utilities::T2pc(T2, CRAMS::H1);
  CHECK(approx(ratio, pc1 / pc2, 1e-9));
}

void test_adiabatic_exact_formula() {
  // dEdX_adiabatic = -(2*v_A)/(3*mu*c) * pc — verify the constructor stores params correctly
  CRAMS::Input in = makeInput();
  CRAMS::Losses L(CRAMS::H1, in);
  const double T = 1. * CRAMS::CGS::GeV;
  const double factorAdv = 2. * in.v_A() / 3. / in.mu() / CRAMS::CGS::cLight;
  const double pc = CRAMS::Utilities::T2pc(T, CRAMS::H1);
  CHECK(approx(L.dEdX_adiabatic(T), -factorAdv * pc, 1e-9));
}

void test_adiabatic_scales_with_vA() {
  // m_factorAdv ∝ v_A → doubling v_A doubles dEdX_adiabatic
  CRAMS::Input in1 = makeInput();
  CRAMS::Input in2 = makeInput();
  const double vA_default_kms = in1.v_A() / (CRAMS::CGS::km / CRAMS::CGS::sec);
  in2.setParam("vA", 2. * vA_default_kms);
  CRAMS::Losses L1(CRAMS::H1, in1);
  CRAMS::Losses L2(CRAMS::H1, in2);
  const double T = 10. * CRAMS::CGS::GeV;
  CHECK(approx(L2.dEdX_adiabatic(T), 2. * L1.dEdX_adiabatic(T), 1e-9));
}

// --- dEdX_ionization: Z² / A scaling ---

void test_ionization_Z2_over_A_scaling() {
  // dEdX_ion ∝ Z²/(A*β²) × betheBlochLog(T)
  // For same β (same T/nucleon), the ratio between two particles is (Z1²/A1)/(Z2²/A2)
  // Use proton (Z=1,A=1) and carbon (Z=6,A=12): ratio = (1/1)/(36/12) = 1/3
  CRAMS::Input in = makeInput();
  CRAMS::Losses L_p(CRAMS::H1, in);
  CRAMS::Losses L_C(CRAMS::C12, in);
  const double T = 10. * CRAMS::CGS::GeV;
  const double ratio = L_p.dEdX_ionization(T) / L_C.dEdX_ionization(T);
  // (Z_p²/A_p) / (Z_C²/A_C) = (1/1) / (36/12) = 1/3
  const double expected = (1. * 1. / 1.) / (6. * 6. / 12.);
  // betheBlochLog has a small A-dependent Q_max correction (~0.04%); relax tolerance
  CHECK(approx(ratio, expected, 1e-3));
}

void test_ionization_increases_with_Z() {
  // Higher Z → more ionization loss at same speed
  CRAMS::Input in = makeInput();
  CRAMS::Losses L_H(CRAMS::H1, in);
  CRAMS::Losses L_C(CRAMS::C12, in);
  const double T = 1. * CRAMS::CGS::GeV;
  CHECK(std::abs(L_C.dEdX_ionization(T)) > std::abs(L_H.dEdX_ionization(T)));
}

void test_ionization_shows_relativistic_rise() {
  // At very high energies: ionization losses increase (relativistic rise ∝ ln γ)
  CRAMS::Losses L(CRAMS::H1, makeInput());
  CHECK(std::abs(L.dEdX_ionization(1. * CRAMS::CGS::TeV)) > std::abs(L.dEdX_ionization(10. * CRAMS::CGS::GeV)));
}

// --- dTdt_ionization ---

void test_dTdt_ionization_linear_in_nH() {
  CRAMS::Losses L(CRAMS::H1, makeInput());
  const double T = 1. * CRAMS::CGS::GeV;
  const double n1 = 0.5 / CRAMS::CGS::cm3;
  const double n2 = 2.0 / CRAMS::CGS::cm3;
  CHECK(approx(L.dTdt_ionization(T, n2), 4. * L.dTdt_ionization(T, n1), 1e-9));
}

void test_dTdt_ionization_Z2_over_A_scaling() {
  CRAMS::Input in = makeInput();
  CRAMS::Losses L_p(CRAMS::H1, in);
  CRAMS::Losses L_C(CRAMS::C12, in);
  const double T = 5. * CRAMS::CGS::GeV;
  const double n_H = 1. / CRAMS::CGS::cm3;
  const double ratio = L_p.dTdt_ionization(T, n_H) / L_C.dTdt_ionization(T, n_H);
  const double expected = (1. * 1. / 1.) / (6. * 6. / 12.);
  CHECK(approx(ratio, expected, 1e-3));
}

// --- getDerivative vs. finite difference ---

void test_getDerivative_matches_finite_difference() {
  // Compare GSL central-difference derivative against manual finite difference
  CRAMS::Losses L(CRAMS::H1, makeInput());
  const double T = 10. * CRAMS::CGS::GeV;
  const double h = 1e-4 * T;
  const double fd = (L.get(T + h) - L.get(T - h)) / (2. * h);
  CHECK(approx(L.getDerivative(T), fd, 1e-4));
}

void test_getDerivative_at_high_energy() {
  CRAMS::Losses L(CRAMS::H1, makeInput());
  const double T = 1. * CRAMS::CGS::TeV;
  const double h = 1e-4 * T;
  const double fd = (L.get(T + h) - L.get(T - h)) / (2. * h);
  CHECK(approx(L.getDerivative(T), fd, 1e-4));
}

// --- adiabatic vs ionization dominance ---

void test_adiabatic_dominates_at_very_high_energy() {
  // At ultra-high energy: dEdX_adiabatic ∝ T, ionization ∝ ln(T) → adiabatic wins
  CRAMS::Losses L(CRAMS::H1, makeInput());
  const double T = 100. * CRAMS::CGS::TeV;
  CHECK(std::abs(L.dEdX_adiabatic(T)) > std::abs(L.dEdX_ionization(T)));
}

void test_ionization_dominates_at_low_energy() {
  // At low energy: ionization losses ∝ 1/β² diverge while adiabatic ∝ pc → 0
  CRAMS::Losses L(CRAMS::H1, makeInput());
  const double T = 0.1 * CRAMS::CGS::GeV;
  CHECK(std::abs(L.dEdX_ionization(T)) > std::abs(L.dEdX_adiabatic(T)));
}

int main() {
  test_get_is_negative();
  test_dEdX_adiabatic_is_negative();
  test_dEdX_ionization_is_negative();
  test_dTdt_ionization_is_positive();
  test_get_equals_sum_of_components();

  test_adiabatic_ratio_equals_momentum_ratio();
  test_adiabatic_exact_formula();
  test_adiabatic_scales_with_vA();

  test_ionization_Z2_over_A_scaling();
  test_ionization_increases_with_Z();
  test_ionization_shows_relativistic_rise();

  test_dTdt_ionization_linear_in_nH();
  test_dTdt_ionization_Z2_over_A_scaling();

  test_getDerivative_matches_finite_difference();
  test_getDerivative_at_high_energy();

  test_adiabatic_dominates_at_very_high_energy();
  test_ionization_dominates_at_low_energy();

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
