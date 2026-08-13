#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "crams/core/cgs.h"
#include "crams/core/pid.h"
#include "crams/secondary.h"

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

static bool approx(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol * std::abs(b) + tol; }

// Power-law grid: Q_i = C * T_i^(-alpha)
static std::vector<double> powerLawQ(const std::vector<double>& T, double C, double alpha) {
  std::vector<double> Q(T.size());
  for (size_t i = 0; i < T.size(); ++i) Q[i] = C * std::pow(T[i], -alpha);
  return Q;
}

// Logarithmically spaced grid
static std::vector<double> logGrid(double Tmin, double Tmax, size_t N) {
  std::vector<double> v(N);
  const double logMin = std::log(Tmin);
  const double logMax = std::log(Tmax);
  for (size_t i = 0; i < N; ++i) v[i] = std::exp(logMin + static_cast<double>(i) * (logMax - logMin) / (N - 1));
  return v;
}

// Shared test grid
static const std::vector<double> T3 = {1. * CRAMS::CGS::GeV, 10. * CRAMS::CGS::GeV, 100. * CRAMS::CGS::GeV};
static const double alpha = 2.7;
static const double C = 1.5e-3;
static const std::vector<double> Q3 = powerLawQ(T3, C, alpha);

// --- boundary behaviour ---

void test_get_zero_below_range() {
  CRAMS::SecondarySource S(CRAMS::H1, T3, Q3);
  CHECK(S.get(0.5 * CRAMS::CGS::GeV) == 0.);
  CHECK(S.get(0.1 * CRAMS::CGS::GeV) == 0.);
}

void test_get_zero_above_range() {
  CRAMS::SecondarySource S(CRAMS::H1, T3, Q3);
  CHECK(S.get(200. * CRAMS::CGS::GeV) == 0.);
  CHECK(S.get(1. * CRAMS::CGS::TeV) == 0.);
}

void test_get_zero_at_lower_boundary() {
  CRAMS::SecondarySource S(CRAMS::H1, T3, Q3);
  CHECK(S.get(T3.front()) == 0.);
}

void test_get_zero_at_upper_boundary() {
  CRAMS::SecondarySource S(CRAMS::H1, T3, Q3);
  CHECK(S.get(T3.back()) == 0.);
}

// --- interpolation correctness ---

void test_get_positive_in_range() {
  CRAMS::SecondarySource S(CRAMS::H1, T3, Q3);
  CHECK(S.get(3. * CRAMS::CGS::GeV) > 0.);
  CHECK(S.get(50. * CRAMS::CGS::GeV) > 0.);
}

void test_get_power_law_exact() {
  // LinearInterpolatorLog is exact for power laws (linearly interpolates in log-log).
  // At any T in the interior, the result must match C * T^(-alpha) to machine precision.
  const auto T_large = logGrid(1. * CRAMS::CGS::GeV, 100. * CRAMS::CGS::GeV, 20);
  const auto Q_large = powerLawQ(T_large, C, alpha);
  CRAMS::SecondarySource S(CRAMS::H1, T_large, Q_large);

  // Test at geometric mid-points between consecutive knots
  for (size_t i = 0; i + 1 < T_large.size(); ++i) {
    const double T_mid = std::sqrt(T_large[i] * T_large[i + 1]);
    const double expected = C * std::pow(T_mid, -alpha);
    CHECK(approx(S.get(T_mid), expected, 1e-12));
  }
}

void test_log_lerp_between_zeros() {
  const auto T = logGrid(1. * CRAMS::CGS::GeV, 100. * CRAMS::CGS::GeV, 20);
  std::vector<double> Q;
  std::fill(Q.begin(), Q.end(), 0.0);
  CRAMS::SecondarySource S(CRAMS::H1, T, Q);

  CHECK(approx(S.get(T[3]), 0.0, 1e-12));                    // on grid
  CHECK(approx(S.get(1.03 * CRAMS::CGS::GeV), 0.0, 1e-12));  // off grid
}

void test_get_matches_exact_at_near_knot() {
  // Spot-check a simple 3-point power-law grid at the geometric midpoint of [1,10] GeV
  CRAMS::SecondarySource S(CRAMS::H1, T3, Q3);
  const double T_mid = std::sqrt(T3[0] * T3[1]);  // sqrt(1*10) = sqrt(10) GeV
  const double expected = C * std::pow(T_mid, -alpha);
  CHECK(approx(S.get(T_mid), expected, 1e-12));
}

void test_get_decreases_with_energy() {
  // Q ∝ T^(-2.7): larger T → smaller source
  CRAMS::SecondarySource S(CRAMS::H1, T3, Q3);
  CHECK(S.get(5. * CRAMS::CGS::GeV) > S.get(20. * CRAMS::CGS::GeV));
}

void test_get_ratio_matches_power_law() {
  // get(T1)/get(T2) = (T1/T2)^(-alpha) for a power-law grid
  CRAMS::SecondarySource S(CRAMS::H1, T3, Q3);
  const double T1 = 2. * CRAMS::CGS::GeV;
  const double T2 = 8. * CRAMS::CGS::GeV;
  const double ratio = S.get(T1) / S.get(T2);
  const double expected = std::pow(T1 / T2, -alpha);
  CHECK(approx(ratio, expected, 1e-12));
}

// --- constructor validation ---

void test_constructor_throws_for_negative_Q() {
  std::vector<double> Q_bad = Q3;
  Q_bad[1] = -1.;
  CHECK_THROW(CRAMS::SecondarySource(CRAMS::H1, T3, Q_bad), std::runtime_error);
}

void test_constructor_throws_for_nan_Q() {
  std::vector<double> Q_bad = Q3;
  Q_bad[0] = std::numeric_limits<double>::quiet_NaN();
  CHECK_THROW(CRAMS::SecondarySource(CRAMS::H1, T3, Q_bad), std::runtime_error);
}

void test_constructor_throws_for_inf_Q() {
  std::vector<double> Q_bad = Q3;
  Q_bad[2] = std::numeric_limits<double>::infinity();
  CHECK_THROW(CRAMS::SecondarySource(CRAMS::H1, T3, Q_bad), std::runtime_error);
}

void test_constructor_accepts_zero_Q_entry() {
  // Zero is valid (no source contribution at that energy)
  std::vector<double> Q_zero = Q3;
  Q_zero[1] = 0.;
  CRAMS::SecondarySource S(CRAMS::H1, T3, Q_zero);
  CHECK(S.get(0.5 * CRAMS::CGS::GeV) == 0.);  // outside range still 0
}

// --- different PIDs ---

void test_get_works_for_heavy_nucleus() {
  CRAMS::SecondarySource S(CRAMS::C12, T3, Q3);
  CHECK(S.get(3. * CRAMS::CGS::GeV) > 0.);
}

void test_get_pid_does_not_affect_interpolation() {
  // PID is stored for logging only; it does not change the interpolated value
  CRAMS::SecondarySource S_p(CRAMS::H1, T3, Q3);
  CRAMS::SecondarySource S_C(CRAMS::C12, T3, Q3);
  const double T = 5. * CRAMS::CGS::GeV;
  CHECK(approx(S_p.get(T), S_C.get(T)));
}

int main() {
  test_get_zero_below_range();
  test_get_zero_above_range();
  test_get_zero_at_lower_boundary();
  test_get_zero_at_upper_boundary();

  test_get_positive_in_range();
  test_get_power_law_exact();
  test_get_matches_exact_at_near_knot();
  test_get_decreases_with_energy();
  test_get_ratio_matches_power_law();

  test_constructor_throws_for_negative_Q();
  test_constructor_throws_for_nan_Q();
  test_constructor_throws_for_inf_Q();
  test_constructor_accepts_zero_Q_entry();

  test_get_works_for_heavy_nucleus();
  test_get_pid_does_not_affect_interpolation();

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
