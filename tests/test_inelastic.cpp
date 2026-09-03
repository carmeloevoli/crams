#include <cmath>
#include <functional>
#include <iostream>

#include "crams/core/cgs.h"
#include "crams/core/pid.h"
#include "crams/inelastic.h"

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

// --- sigma_pp ---

void test_sigma_pp_is_zero_below_threshold() {
  CHECK(CRAMS::sigma_pp(0.1 * CRAMS::CGS::GeV) == 0.);
  CHECK(CRAMS::sigma_pp(0.27 * CRAMS::CGS::GeV) == 0.);
}

void test_sigma_pp_is_positive_above_threshold() {
  CHECK(CRAMS::sigma_pp(1. * CRAMS::CGS::GeV) > 0.);
  CHECK(CRAMS::sigma_pp(100. * CRAMS::CGS::GeV) > 0.);
}

void test_sigma_pp_physically_reasonable_at_10GeV() {
  // pp inelastic cross-section near 10 GeV is ~30 mbarn from the Tan & Ng parameterisation
  const double sigma = CRAMS::sigma_pp(10. * CRAMS::CGS::GeV);
  CHECK(sigma > 20. * CRAMS::CGS::mbarn && sigma < 45. * CRAMS::CGS::mbarn);
}

void test_sigma_pp_rises_above_threshold() {
  const double s_low = CRAMS::sigma_pp(0.35 * CRAMS::CGS::GeV);
  const double s_high = CRAMS::sigma_pp(1.0 * CRAMS::CGS::GeV);
  CHECK(s_high > s_low);
}

// --- InelasticXsecST98 ---

void test_st98_H1_matches_sigma_pp() {
  CRAMS::InelasticXsecST98 xsec;
  const double T = 10. * CRAMS::CGS::GeV;
  CHECK(approx(xsec.getXsecOnHtarget(CRAMS::H1, T), CRAMS::sigma_pp(T), 1e-9));
}

void test_st98_nuclei_match_sigma_ST() {
  CRAMS::InelasticXsecST98 xsec;
  const double T = 10. * CRAMS::CGS::GeV;
  CHECK(approx(xsec.getXsecOnHtarget(CRAMS::C12, T), CRAMS::sigma_ST(T, CRAMS::C12.getA()), 1e-9));
  CHECK(approx(xsec.getXsecOnHtarget(CRAMS::Fe56, T), CRAMS::sigma_ST(T, CRAMS::Fe56.getA()), 1e-9));
}

void test_st98_ISM_xsec_exact_factor() {
  CRAMS::InelasticXsecST98 xsec;
  const double T = 10. * CRAMS::CGS::GeV;
  const double expected = (1. + CRAMS::CGS::K_He * CRAMS::CGS::f_He) / (1. + CRAMS::CGS::f_He);
  CHECK(approx(xsec.getXsecOnISM(CRAMS::C12, T) / xsec.getXsecOnHtarget(CRAMS::C12, T), expected, 1e-9));
}

// --- InXsecTripathi99: proton (H1) ---

void test_tripathi_H1_matches_sigma_pp() {
  // For the proton, getXsecOnHtarget delegates exactly to sigma_pp
  CRAMS::InXsecTripathi99 xsec;
  const double T = 10. * CRAMS::CGS::GeV;
  CHECK(approx(xsec.getXsecOnHtarget(CRAMS::H1, T), CRAMS::sigma_pp(T), 1e-9));
}

void test_tripathi_H1_below_threshold_is_clamped() {
  // sigma_pp = 0 below threshold; getXsecOnHtarget clamps to 1e-10 mbarn
  CRAMS::InXsecTripathi99 xsec;
  CHECK(xsec.getXsecOnHtarget(CRAMS::H1, 0.1 * CRAMS::CGS::GeV) >= 1e-10 * CRAMS::CGS::mbarn);
}

// --- InXsecTripathi99: nuclei ---

void test_tripathi_He4_positive_in_GCR_range() {
  CRAMS::InXsecTripathi99 xsec;
  CHECK(xsec.getXsecOnHtarget(CRAMS::He4, 1. * CRAMS::CGS::GeV) > 0.);
  CHECK(xsec.getXsecOnHtarget(CRAMS::He4, 100. * CRAMS::CGS::GeV) > 0.);
}

void test_tripathi_sigma_monotonic_in_A() {
  // Heavier projectiles have larger geometric cross-sections: He4 < C12 < Fe56
  CRAMS::InXsecTripathi99 xsec;
  const double T = 10. * CRAMS::CGS::GeV;
  CHECK(xsec.getXsecOnHtarget(CRAMS::C12, T) > xsec.getXsecOnHtarget(CRAMS::He4, T));
  CHECK(xsec.getXsecOnHtarget(CRAMS::Fe56, T) > xsec.getXsecOnHtarget(CRAMS::C12, T));
}

void test_tripathi_He4_value_at_10GeV() {
  // He4 on H at ~10 GeV/n: table gives ~112 mbarn
  CRAMS::InXsecTripathi99 xsec;
  const double sigma = xsec.getXsecOnHtarget(CRAMS::He4, 10. * CRAMS::CGS::GeV);
  CHECK(sigma > 50. * CRAMS::CGS::mbarn && sigma < 200. * CRAMS::CGS::mbarn);
}

void test_tripathi_C12_value_at_10GeV() {
  // C12 on H at ~10 GeV/n: table gives ~251 mbarn
  CRAMS::InXsecTripathi99 xsec;
  const double sigma = xsec.getXsecOnHtarget(CRAMS::C12, 10. * CRAMS::CGS::GeV);
  CHECK(sigma > 100. * CRAMS::CGS::mbarn && sigma < 500. * CRAMS::CGS::mbarn);
}

void test_tripathi_Fe56_value_at_10GeV() {
  // Fe56 on H at ~10 GeV/n: table gives ~719 mbarn
  CRAMS::InXsecTripathi99 xsec;
  const double sigma = xsec.getXsecOnHtarget(CRAMS::Fe56, 10. * CRAMS::CGS::GeV);
  CHECK(sigma > 300. * CRAMS::CGS::mbarn && sigma < 1500. * CRAMS::CGS::mbarn);
}

static bool throwsRuntimeError(const std::function<void()>& fn) {
  try {
    fn();
  } catch (const std::runtime_error&) {
    return true;
  } catch (...) {
    return false;
  }
  return false;
}

void test_tripathi_high_energy_extrapolates_out_of_range() {
  // T above the table maximum (1e5 GeV) is out of range and must throw
  CRAMS::InXsecTripathi99 xsec;
  CHECK(!throwsRuntimeError([&] { xsec.getXsecOnHtarget(CRAMS::C12, 1e6 * CRAMS::CGS::GeV); }));
}

void test_tripathi_low_energy_throws_out_of_range() {
  // T below the table minimum (0.01 GeV) is out of range and must throw
  CRAMS::InXsecTripathi99 xsec;
  CHECK(throwsRuntimeError([&] { xsec.getXsecOnHtarget(CRAMS::C12, 1e-3 * CRAMS::CGS::GeV); }));
}

void test_tripathi_at_table_bounds_does_not_throw() {
  // Exactly at m_T_min (0.01 GeV) and m_T_max (1e5 GeV) is in range and must succeed
  CRAMS::InXsecTripathi99 xsec;
  CHECK(!throwsRuntimeError([&] { xsec.getXsecOnHtarget(CRAMS::C12, 0.01 * CRAMS::CGS::GeV); }));
  CHECK(!throwsRuntimeError([&] { xsec.getXsecOnHtarget(CRAMS::C12, 1e5 * CRAMS::CGS::GeV); }));
}

// --- InXsecGlauber ---

void test_glauber_loads_and_returns_positive() {
  // Glauber shares the table machinery; constructing it must load the file and
  // return positive cross-sections for nuclei in the GCR range.
  CRAMS::InXsecGlauber xsec;
  CHECK(xsec.getXsecOnHtarget(CRAMS::C12, 10. * CRAMS::CGS::GeV) > 0.);
  CHECK(xsec.getXsecOnHtarget(CRAMS::Fe56, 10. * CRAMS::CGS::GeV) > 0.);
}

void test_glauber_sigma_monotonic_in_A() {
  // Heavier projectiles have larger inelastic cross-sections
  CRAMS::InXsecGlauber xsec;
  const double T = 10. * CRAMS::CGS::GeV;
  CHECK(xsec.getXsecOnHtarget(CRAMS::Fe56, T) > xsec.getXsecOnHtarget(CRAMS::C12, T));
}

void test_glauber_low_energy_throws_out_of_range() {
  CRAMS::InXsecGlauber xsec;
  CHECK(!throwsRuntimeError([&] { xsec.getXsecOnHtarget(CRAMS::C12, 1e6 * CRAMS::CGS::GeV); }));
  CHECK(throwsRuntimeError([&] { xsec.getXsecOnHtarget(CRAMS::C12, 1e-3 * CRAMS::CGS::GeV); }));
}

// --- InelasticXsec::getXsecOnISM ---

void test_tripathi_ISM_xsec_exact_factor() {
  // getXsecOnISM = sigma_H * (1 + K_He * f_He) / (1 + f_He)
  CRAMS::InXsecTripathi99 xsec;
  const double T = 10. * CRAMS::CGS::GeV;
  const double expected = (1. + CRAMS::CGS::K_He * CRAMS::CGS::f_He) / (1. + CRAMS::CGS::f_He);
  CHECK(approx(xsec.getXsecOnISM(CRAMS::C12, T) / xsec.getXsecOnHtarget(CRAMS::C12, T), expected, 1e-9));
}

void test_tripathi_ISM_greater_than_H_target() {
  // K_He ≈ 2.52 > 1 and f_He = 0.08 > 0, so ISM factor > 1
  CRAMS::InXsecTripathi99 xsec;
  const double T = 10. * CRAMS::CGS::GeV;
  CHECK(xsec.getXsecOnISM(CRAMS::He4, T) > xsec.getXsecOnHtarget(CRAMS::He4, T));
}

int main() {
  test_sigma_pp_is_zero_below_threshold();
  test_sigma_pp_is_positive_above_threshold();
  test_sigma_pp_physically_reasonable_at_10GeV();
  test_sigma_pp_rises_above_threshold();

  test_st98_H1_matches_sigma_pp();
  test_st98_nuclei_match_sigma_ST();
  test_st98_ISM_xsec_exact_factor();

  test_tripathi_H1_matches_sigma_pp();
  test_tripathi_H1_below_threshold_is_clamped();
  test_tripathi_He4_positive_in_GCR_range();
  test_tripathi_sigma_monotonic_in_A();
  test_tripathi_He4_value_at_10GeV();
  test_tripathi_C12_value_at_10GeV();
  test_tripathi_Fe56_value_at_10GeV();
  test_tripathi_high_energy_extrapolates_out_of_range();
  test_tripathi_low_energy_throws_out_of_range();
  test_tripathi_at_table_bounds_does_not_throw();

  test_glauber_loads_and_returns_positive();
  test_glauber_sigma_monotonic_in_A();
  test_glauber_low_energy_throws_out_of_range();

  test_tripathi_ISM_xsec_exact_factor();
  test_tripathi_ISM_greater_than_H_target();

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
