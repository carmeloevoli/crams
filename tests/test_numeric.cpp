#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "crams/utils/numeric.h"

using namespace CRAMS::Numeric;

namespace {

int failures = 0;

void check(bool ok, const char* msg) {
  if (ok)
    std::cout << "PASS: " << msg << "\n";
  else {
    std::cerr << "FAIL: " << msg << "\n";
    ++failures;
  }
}

bool approx(double a, double b, double tol = 1e-6) { return std::abs(a / b - 1.0) < tol; }

// ---------------------------------------------------------------------------
// getLowerIndex
// ---------------------------------------------------------------------------

void test_getLowerIndex() {
  const std::vector<double> v = {1.0, 2.0, 3.0, 4.0, 5.0};

  // Interior intervals
  check(getLowerIndex(v, 1.5) == 0, "getLowerIndex: 1.5 in [1,2] → 0");
  check(getLowerIndex(v, 2.5) == 1, "getLowerIndex: 2.5 in [2,3] → 1");
  check(getLowerIndex(v, 4.5) == 3, "getLowerIndex: 4.5 in [4,5] → 3");

  // Exact grid points
  check(getLowerIndex(v, 1.0) == 0, "getLowerIndex: front → 0");
  check(getLowerIndex(v, 2.0) == 1, "getLowerIndex: 2.0 (exact grid) → 1");
  check(getLowerIndex(v, 3.0) == 2, "getLowerIndex: 3.0 (exact grid) → 2");
  check(getLowerIndex(v, 5.0) == 3, "getLowerIndex: back → size-2");

  // Result always gives valid i and i+1
  for (double x : {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0}) {
    const size_t i = getLowerIndex(v, x);
    check(i + 1 < v.size(), "getLowerIndex: i+1 always valid");
    check(v[i] <= x && x <= v[i + 1], "getLowerIndex: v[i] <= x <= v[i+1]");
  }
}

// ---------------------------------------------------------------------------
// LinearInterpolator
// ---------------------------------------------------------------------------

void test_LinearInterpolator() {
  const std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
  const std::vector<double> y = {0.0, 1.0, 2.0, 3.0, 4.0};  // identity: y = x

  // Exact at grid points
  check(approx(LinearInterpolator(x, y, 0.0), 0.0, 1e-12) || LinearInterpolator(x, y, 0.0) == 0.0,
        "LinearInterp: front grid point");
  check(approx(LinearInterpolator(x, y, 4.0), 4.0), "LinearInterp: back grid point");
  check(approx(LinearInterpolator(x, y, 2.0), 2.0), "LinearInterp: mid grid point");

  // Previously buggy: x_new in first interval [x[0], x[1]]
  check(approx(LinearInterpolator(x, y, 0.5), 0.5), "LinearInterp: first interval (was buggy)");
  check(approx(LinearInterpolator(x, y, 0.25), 0.25), "LinearInterp: first interval 0.25");

  // Interior midpoints
  check(approx(LinearInterpolator(x, y, 1.5), 1.5), "LinearInterp: midpoint 1.5");
  check(approx(LinearInterpolator(x, y, 3.7), 3.7), "LinearInterp: 3.7");

  // Non-trivial linear function: y = 2x + 1
  const std::vector<double> y2 = {1.0, 3.0, 5.0, 7.0, 9.0};
  check(approx(LinearInterpolator(x, y2, 1.5), 4.0), "LinearInterp: y=2x+1 at 1.5");
  check(approx(LinearInterpolator(x, y2, 0.3), 1.6), "LinearInterp: y=2x+1 at 0.3");

  // Out-of-range throws
  try {
    LinearInterpolator(x, y, -0.1);
    check(false, "LinearInterp: below range should throw");
  } catch (const std::invalid_argument&) {
    check(true, "LinearInterp: below range throws");
  }
  try {
    LinearInterpolator(x, y, 4.1);
    check(false, "LinearInterp: above range should throw");
  } catch (const std::invalid_argument&) {
    check(true, "LinearInterp: above range throws");
  }
}

// ---------------------------------------------------------------------------
// LinearInterpolatorLog
// ---------------------------------------------------------------------------

void test_LinearInterpolatorLog() {
  // Power law y = x^2 is exact in log-log space
  const std::vector<double> x = {1.0, 10.0, 100.0, 1000.0};
  const std::vector<double> y = {1.0, 100.0, 10000.0, 1000000.0};

  // Exact at grid points
  check(approx(LinearInterpolatorLog(x, y, 1.0), 1.0), "LinearInterpLog: x=1 grid point");
  check(approx(LinearInterpolatorLog(x, y, 10.0), 100.0), "LinearInterpLog: x=10 grid point");
  check(approx(LinearInterpolatorLog(x, y, 1000.0), 1e6), "LinearInterpLog: x=1000 back");

  // Midpoint in log space: x=sqrt(10), y=10
  check(approx(LinearInterpolatorLog(x, y, std::sqrt(10.0)), 10.0, 1e-5),
        "LinearInterpLog: x=sqrt(10) → y=10 (exact power law)");

  // Another midpoint: x=sqrt(1000) ≈ 31.62, y=1000
  check(approx(LinearInterpolatorLog(x, y, std::sqrt(1000.0)), 1000.0, 1e-5),
        "LinearInterpLog: x=sqrt(1000) → y=1000");

  // First interval (was also potentially buggy with old getLowerIndex)
  const double x_first = std::sqrt(1.0 * 10.0);  // geometric mean of first two points
  check(approx(LinearInterpolatorLog(x, y, x_first), x_first * x_first, 1e-5),
        "LinearInterpLog: first interval midpoint");

  // Out-of-range throws
  try {
    LinearInterpolatorLog(x, y, 0.5);
    check(false, "LinearInterpLog: below range should throw");
  } catch (const std::invalid_argument&) {
    check(true, "LinearInterpLog: below range throws");
  }
}

// ---------------------------------------------------------------------------
// QAGIntegration
// ---------------------------------------------------------------------------

void test_QAGIntegration() {
  // int_0^1 x^2 dx = 1/3
  auto poly = [](double x) { return x * x; };
  check(approx(QAGIntegration<double>(poly, 0.0, 1.0), 1.0 / 3.0),
        "QAG: integral of x^2 on [0,1] = 1/3");

  // int_0^pi sin(x) dx = 2
  auto sinf = [](double x) { return std::sin(x); };
  check(approx(QAGIntegration<double>(sinf, 0.0, M_PI), 2.0),
        "QAG: integral of sin on [0,pi] = 2");

  // int_1^e 1/x dx = 1 (ln(e) - ln(1) = 1)
  auto invx = [](double x) { return 1.0 / x; };
  check(approx(QAGIntegration<double>(invx, 1.0, M_E), 1.0),
        "QAG: integral of 1/x on [1,e] = 1");
}

// ---------------------------------------------------------------------------
// QAGSIntegration
// ---------------------------------------------------------------------------

void test_QAGSIntegration() {
  auto poly = [](double x) { return x * x; };
  check(approx(QAGSIntegration<double>(poly, 0.0, 1.0), 1.0 / 3.0),
        "QAGS: integral of x^2 on [0,1] = 1/3");

  // QAGS handles integrable singularities: int_0^1 1/sqrt(x) dx = 2
  auto sqrtInv = [](double x) { return 1.0 / std::sqrt(x + 1e-10); };
  check(approx(QAGSIntegration<double>(sqrtInv, 0.0, 1.0), 2.0, 1e-3),
        "QAGS: integral of 1/sqrt(x) on [0,1] ≈ 2");
}

// ---------------------------------------------------------------------------
// simpsonIntegration
// ---------------------------------------------------------------------------

void test_simpsonIntegration() {
  // int_0^1 x^2 dx = 1/3 (Simpson's is exact for polynomials deg <= 3)
  auto poly = [](double x) { return x * x; };
  check(approx(simpsonIntegration<double>(poly, 0.0, 1.0, 100), 1.0 / 3.0),
        "Simpson: x^2 on [0,1]");

  // int_0^1 x^3 dx = 1/4 (exact for Simpson with even N)
  auto cubic = [](double x) { return x * x * x; };
  check(approx(simpsonIntegration<double>(cubic, 0.0, 1.0, 100), 1.0 / 4.0),
        "Simpson: x^3 on [0,1]");

  // int_0^pi sin(x) dx = 2 (converges with enough points)
  auto sinf = [](double x) { return std::sin(x); };
  check(approx(simpsonIntegration<double>(sinf, 0.0, M_PI, 1000), 2.0, 1e-5),
        "Simpson: sin on [0,pi]");

  // Odd N is auto-corrected to even
  check(approx(simpsonIntegration<double>(poly, 0.0, 1.0, 99), 1.0 / 3.0),
        "Simpson: odd N auto-corrected");
}

// ---------------------------------------------------------------------------
// interpolate2d
// ---------------------------------------------------------------------------

void test_interpolate2d() {
  // f(x,y) = x + y — bilinear interpolation is exact for this
  const std::vector<double> x = {0.0, 1.0, 2.0};
  const std::vector<double> y = {0.0, 1.0, 2.0};
  // z[j + ny*i] = f(x[i], y[j]) = x[i] + y[j]
  std::vector<double> z(9);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      z[j + 3 * i] = x[i] + y[j];

  check(approx(interpolate2d(x, y, z, 0.5, 0.5), 1.0), "interpolate2d: f=x+y at (0.5,0.5)=1");
  check(approx(interpolate2d(x, y, z, 1.0, 1.0), 2.0), "interpolate2d: f=x+y at (1,1)=2");
  check(std::abs(interpolate2d(x, y, z, 0.0, 0.0)) < 1e-10, "interpolate2d: corner (0,0)=0");
  check(approx(interpolate2d(x, y, z, 2.0, 2.0), 4.0), "interpolate2d: corner (2,2)=4");
  check(approx(interpolate2d(x, y, z, 1.5, 0.3), 1.8), "interpolate2d: f=x+y at (1.5,0.3)=1.8");
}

}  // namespace

int main() {
  test_getLowerIndex();
  test_LinearInterpolator();
  test_LinearInterpolatorLog();
  test_QAGIntegration();
  test_QAGSIntegration();
  test_simpsonIntegration();
  test_interpolate2d();

  if (failures == 0)
    std::cout << "\nAll " << __FILE__ << " tests passed.\n";
  else
    std::cerr << "\n" << failures << " test(s) failed.\n";

  return failures > 0 ? 1 : 0;
}
