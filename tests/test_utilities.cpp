#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "crams/core/cgs.h"
#include "crams/utils/utilities.h"

using namespace CRAMS;
using namespace CRAMS::Utilities;

// ---------------------------------------------------------------------------
// Compile-time checks
// ---------------------------------------------------------------------------

static_assert(pow2(3) == 9, "pow2 int");
static_assert(pow3(2) == 8, "pow3 int");
static_assert(pow4(2) == 16, "pow4 int");
static_assert(pow2(3.0) == 9.0, "pow2 double");
static_assert(pow2(-4) == 16, "pow2 negative");
static_assert(pow4(3) == pow2(pow2(3)), "pow4 = pow2 of pow2");

// ---------------------------------------------------------------------------
// Runtime helpers
// ---------------------------------------------------------------------------

namespace {

int failures = 0;

void check(bool ok, const char* msg) {
  if (ok) {
    std::cout << "PASS: " << msg << "\n";
  } else {
    std::cerr << "FAIL: " << msg << "\n";
    ++failures;
  }
}

bool approx(double a, double b, double tol = 1e-5) { return std::abs(a / b - 1.0) < tol; }

// ---------------------------------------------------------------------------
// T2beta
// ---------------------------------------------------------------------------

void test_T2beta() {
  check(T2beta(0.0) == 0.0, "T2beta(0) = 0");
  // At T = mp*c^2: beta = sqrt(3)/2
  check(approx(T2beta(CGS::protonMassC2), std::sqrt(3.) / 2.), "T2beta(mp) = sqrt(3)/2");
  // Large T: beta approaches 1
  check(T2beta(1e6 * CGS::protonMassC2) < 1.0, "T2beta < 1 always");
  check(T2beta(1e6 * CGS::protonMassC2) > 0.9999, "T2beta(large T) -> 1");
  // Throws on negative T
  try {
    T2beta(-1.0);
    check(false, "T2beta(-1) should throw");
  } catch (const std::invalid_argument&) {
    check(true, "T2beta(-1) throws invalid_argument");
  }
}

// ---------------------------------------------------------------------------
// T2gamma
// ---------------------------------------------------------------------------

void test_T2gamma() {
  check(T2gamma(0.0) == 1.0, "T2gamma(0) = 1");
  // At T = mp*c^2: gamma = 2
  check(approx(T2gamma(CGS::protonMassC2), 2.0), "T2gamma(mp) = 2");
  // gamma always >= 1
  check(T2gamma(CGS::GeV) >= 1.0, "T2gamma >= 1");
  // beta and gamma satisfy beta = sqrt(1 - 1/gamma^2)
  const double T = 10. * CGS::GeV;
  const double beta = T2beta(T);
  const double gamma = T2gamma(T);
  check(approx(beta, std::sqrt(1. - 1. / pow2(gamma))), "beta = sqrt(1 - 1/gamma^2)");
  // Throws on negative T
  try {
    T2gamma(-1.0);
    check(false, "T2gamma(-1) should throw");
  } catch (const std::invalid_argument&) {
    check(true, "T2gamma(-1) throws invalid_argument");
  }
}

// ---------------------------------------------------------------------------
// T2pc and R2T
// ---------------------------------------------------------------------------

void test_T2pc_R2T() {
  const PID proton(1, 1);
  const PID carbon(6, 12);

  // T2pc at T=0 gives 0
  check(T2pc(0.0, proton) == 0.0, "T2pc(0, proton) = 0");
  // T2pc grows with T
  check(T2pc(10. * CGS::GeV, proton) > T2pc(1. * CGS::GeV, proton), "T2pc increases with T");
  // Heavier nucleus has larger pc at same T per nucleon
  check(T2pc(1. * CGS::GeV, carbon) > T2pc(1. * CGS::GeV, proton), "T2pc(C) > T2pc(p) at same T");

  // R2T / T2pc roundtrip: R = T2pc(T, pid) / Z  =>  R2T(R, pid) = T
  const double T0 = 10. * CGS::GeV;
  const double R_proton = T2pc(T0, proton) / proton.getZ();
  check(approx(R2T(R_proton, proton), T0), "R2T(T2pc/Z, proton) roundtrip");
  const double R_carbon = T2pc(T0, carbon) / carbon.getZ();
  check(approx(R2T(R_carbon, carbon), T0), "R2T(T2pc/Z, carbon) roundtrip");

  // R2T throws on negative R
  try {
    R2T(-1.0, proton);
    check(false, "R2T(-1) should throw");
  } catch (const std::invalid_argument&) {
    check(true, "R2T(-1) throws invalid_argument");
  }
}

// ---------------------------------------------------------------------------
// LinAxis
// ---------------------------------------------------------------------------

void test_LinAxis() {
  const auto v = LinAxis(0.0, 10.0, 11);
  check(v.size() == 11, "LinAxis size");
  check(v.front() == 0.0, "LinAxis first = min");
  check(v.back() == 10.0, "LinAxis last = max");
  // Equal spacing
  const double dx = v[1] - v[0];
  bool uniform = true;
  for (size_t i = 1; i < v.size(); ++i)
    if (!approx(v[i] - v[i - 1], dx)) uniform = false;
  check(uniform, "LinAxis uniform spacing");
  // Monotonically increasing
  bool mono = true;
  for (size_t i = 1; i < v.size(); ++i)
    if (v[i] <= v[i - 1]) mono = false;
  check(mono, "LinAxis monotonically increasing");
  // Throws on bad args
  try {
    LinAxis(5., 1., 10);
    check(false, "LinAxis min>max should throw");
  } catch (const std::invalid_argument&) {
    check(true, "LinAxis min>max throws");
  }
  try {
    LinAxis(0., 1., 1);
    check(false, "LinAxis size=1 should throw");
  } catch (const std::invalid_argument&) {
    check(true, "LinAxis size=1 throws");
  }
}

// ---------------------------------------------------------------------------
// LogAxis
// ---------------------------------------------------------------------------

void test_LogAxis() {
  const auto v = LogAxis(1.0, 1000.0, 7);
  check(v.size() == 7, "LogAxis size");
  check(approx(v.front(), 1.0), "LogAxis first = min");
  check(approx(v.back(), 1000.0), "LogAxis last = max");
  // Constant log ratio between consecutive elements
  const double ratio = v[1] / v[0];
  bool log_uniform = true;
  for (size_t i = 1; i < v.size(); ++i)
    if (!approx(v[i] / v[i - 1], ratio)) log_uniform = false;
  check(log_uniform, "LogAxis constant log-ratio");
  // Monotonically increasing
  bool mono = true;
  for (size_t i = 1; i < v.size(); ++i)
    if (v[i] <= v[i - 1]) mono = false;
  check(mono, "LogAxis monotonically increasing");
  // All values positive
  check(isGoodAndPositive(v), "LogAxis all positive");
  // Throws on bad args
  try {
    LogAxis(5., 1., 10);
    check(false, "LogAxis min>max should throw");
  } catch (const std::invalid_argument&) {
    check(true, "LogAxis min>max throws");
  }
}

// ---------------------------------------------------------------------------
// isGoodAndPositive
// ---------------------------------------------------------------------------

void test_isGoodAndPositive() {
  check(isGoodAndPositive({1.0, 2.0, 3.0}), "isGoodAndPositive: all positive");
  check(isGoodAndPositive({0.0, 1.0}), "isGoodAndPositive: zero is ok");
  check(!isGoodAndPositive({1.0, -1.0}), "isGoodAndPositive: negative fails");
  check(!isGoodAndPositive({1.0, std::numeric_limits<double>::quiet_NaN()}), "isGoodAndPositive: NaN fails");
  check(!isGoodAndPositive({1.0, std::numeric_limits<double>::infinity()}), "isGoodAndPositive: +inf fails");
  check(!isGoodAndPositive({1.0, -std::numeric_limits<double>::infinity()}), "isGoodAndPositive: -inf fails");
  check(isGoodAndPositive({}), "isGoodAndPositive: empty vector is ok");
}

// ---------------------------------------------------------------------------
// simplifyKey
// ---------------------------------------------------------------------------

void test_simplifyKey() {
  check(simplifyKey("Hello World") == "helloworld", "simplifyKey: spaces removed");
  check(simplifyKey("diffusion_coefficient") == "diffusioncoefficient", "simplifyKey: underscores removed");
  check(simplifyKey("UPPER_CASE") == "uppercase", "simplifyKey: lowercase");
  check(simplifyKey("") == "", "simplifyKey: empty string");
  check(simplifyKey("already") == "already", "simplifyKey: no-op on clean string");
}

// ---------------------------------------------------------------------------
// fileExists
// ---------------------------------------------------------------------------

void test_fileExists() {
  check(!fileExists("/nonexistent/path/to/file.txt"), "fileExists: missing file = false");
  // Write a temp file and check it exists
  const std::string tmp = "/tmp/test_crams_fileexists.txt";
  {
    std::ofstream f(tmp);
    f << "test";
  }
  check(fileExists(tmp), "fileExists: existing file = true");
  std::remove(tmp.c_str());
}

// ---------------------------------------------------------------------------
// inRange
// ---------------------------------------------------------------------------

void test_inRange() {
  check(inRange(5.0, {0.0, 10.0}), "inRange: interior");
  check(inRange(0.0, {0.0, 10.0}), "inRange: at lower bound");
  check(inRange(10.0, {0.0, 10.0}), "inRange: at upper bound");
  check(!inRange(-1.0, {0.0, 10.0}), "inRange: below range");
  check(!inRange(11.0, {0.0, 10.0}), "inRange: above range");
}

// ---------------------------------------------------------------------------
// loadColumn
// ---------------------------------------------------------------------------

void test_loadColumn() {
  const std::string tmp = "/tmp/test_crams_loadcol.txt";
  {
    std::ofstream f(tmp);
    f << "# header\n";
    f << "1.0 2.0 3.0\n";
    f << "4.0 5.0 6.0\n";
    f << "7.0 8.0 9.0\n";
  }
  const auto col0 = loadColumn(tmp, 0, 1);
  check(col0.size() == 3, "loadColumn: correct row count");
  check(approx(col0[0], 1.0) && approx(col0[1], 4.0) && approx(col0[2], 7.0), "loadColumn: col 0 values");
  const auto col1 = loadColumn(tmp, 1, 1);
  check(approx(col1[0], 2.0) && approx(col1[1], 5.0) && approx(col1[2], 8.0), "loadColumn: col 1 values");
  std::remove(tmp.c_str());
  // Throws on missing file
  try {
    loadColumn("/nonexistent/file.txt", 0, 0);
    check(false, "loadColumn: missing file should throw");
  } catch (const std::runtime_error&) {
    check(true, "loadColumn: missing file throws runtime_error");
  }
}

}  // namespace

int main() {
  test_T2beta();
  test_T2gamma();
  test_T2pc_R2T();
  test_LinAxis();
  test_LogAxis();
  test_isGoodAndPositive();
  test_simplifyKey();
  test_fileExists();
  test_inRange();
  test_loadColumn();

  if (failures == 0)
    std::cout << "\nAll " << __FILE__ << " tests passed.\n";
  else
    std::cerr << "\n" << failures << " test(s) failed.\n";

  return failures > 0 ? 1 : 0;
}
