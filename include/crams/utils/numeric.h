#ifndef INCLUDE_NUMERIC_H
#define INCLUDE_NUMERIC_H

#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline2d.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

namespace CRAMS {
namespace Numeric {

// Returns i such that v[i] <= x <= v[i+1]. Binary search: O(log n).
template <typename T>
size_t getLowerIndex(const std::vector<T>& v, T x) {
  assert(x >= v.front() && x <= v.back());
  // upper_bound gives the first element strictly greater than x
  auto it = std::upper_bound(v.begin(), v.end(), x);
  if (it != v.begin()) --it;
  size_t i = static_cast<size_t>(it - v.begin());
  // clamp so i+1 is always a valid index (handles x == v.back())
  if (i >= v.size() - 1) i = v.size() - 2;
  return i;
}

// Linear interpolation in linear-linear space.
template <typename T>
T LinearInterpolator(const std::vector<T>& x, const std::vector<T>& y, T x_new) {
  if (x_new < x.front() || x_new > x.back()) throw std::invalid_argument("x_new out of range in LinearInterpolator");
  const size_t i = getLowerIndex(x, x_new);
  const T t = (x_new - x[i]) / (x[i + 1] - x[i]);
  return y[i] * (1. - t) + y[i + 1] * t;
}

// Linear interpolation in log-log space. Exact for power laws.
template <typename T>
T LinearInterpolatorLog(const std::vector<T>& x, const std::vector<T>& y, T x_new) {
  if (x_new < x.front() || x_new > x.back()) throw std::invalid_argument("x_new out of range in LinearInterpolatorLog");
  const size_t i = getLowerIndex(x, x_new);
  const double t = (std::log(x_new) - std::log(x[i])) / (std::log(x[i + 1]) - std::log(x[i]));
  return std::exp(std::log(y[i]) * (1. - t) + std::log(y[i + 1]) * t);
}

// GSL adaptive integration (QAG) over a finite interval.
template <typename T>
T QAGIntegration(std::function<T(T)> f, T start, T stop, size_t limit = 1000, double rel_error = 1e-4) {
  gsl_function F;
  F.function = [](double x, void* vf) -> double { return (*static_cast<std::function<double(double)>*>(vf))(x); };
  F.params = &f;
  double result, error;
  gsl_integration_workspace* ws = gsl_integration_workspace_alloc(limit);
  gsl_integration_qag(&F, start, stop, 0., rel_error, limit, GSL_INTEG_GAUSS31, ws, &result, &error);
  gsl_integration_workspace_free(ws);
  return T(result);
}

// GSL adaptive integration (QAGS) with singularity handling.
template <typename T>
T QAGSIntegration(std::function<T(T)> f, T start, T stop, size_t limit = 1000, double rel_error = 1e-4) {
  gsl_function F;
  F.function = [](double x, void* vf) -> double { return (*static_cast<std::function<double(double)>*>(vf))(x); };
  F.params = &f;
  double result, error;
  gsl_integration_workspace* ws = gsl_integration_workspace_alloc(limit);
  gsl_integration_qags(&F, start, stop, 0., rel_error, limit, ws, &result, &error);
  gsl_integration_workspace_free(ws);
  return T(result);
}

// GSL adaptive integration (QAGIU) over [start, +∞).
template <typename T>
T QAGIUIntegration(std::function<T(T)> f, T start, size_t limit = 1000, double rel_error = 1e-4) {
  gsl_function F;
  F.function = [](double x, void* vf) -> double { return (*static_cast<std::function<double(double)>*>(vf))(x); };
  F.params = &f;
  double result, error;
  gsl_integration_workspace* ws = gsl_integration_workspace_alloc(limit);
  gsl_integration_qagiu(&F, start, 0., rel_error, limit, ws, &result, &error);
  gsl_integration_workspace_free(ws);
  return T(result);
}

// Dimensionless spectral integral 4 pi ∫ x^(2-slope) (sqrt(1 + x^2) - 1) dx, slope ∈ (4, 5).
// Used to normalise the primary CR source to the SNR energy budget.
inline double gammaIntegral(double slope) {
  if (!(slope > 4.0 && slope < 5.0)) throw std::invalid_argument("slope must be in (4, 5)");

  double result;
  if (slope < 4.1) {
    // Near slope=4 the integrand decays slowly; use semi-infinite quadrature directly.
    const auto integrand = [slope](double x) -> double {
      return std::pow(x, 2. - slope) * (std::sqrt(x * x + 1.) - 1.);
    };
    result = QAGIUIntegration<double>(integrand, 0.);
  } else {
    // Log substitution x = e^y concentrates the integrand and avoids the slow tail.
    const auto integrand = [slope](double y) -> double {
      const double x = std::exp(y);
      return x * std::pow(x, 2. - slope) * (std::sqrt(x * x + 1.) - 1.);
    };
    result = QAGIntegration<double>(integrand, std::log(1e-5), std::log(1e10));
  }
  return 4. * M_PI * result;
}

// Composite Simpson's rule. N is rounded up to the nearest even number.
template <typename T>
T simpsonIntegration(std::function<T(T)> f, T start, T stop, size_t N = 100) {
  if (N % 2 != 0) ++N;
  const T h = (stop - start) / static_cast<T>(N);
  T result = f(start) + f(stop);
  for (size_t i = 1; i < N; ++i) {
    result += static_cast<T>(i % 2 == 0 ? 2 : 4) * f(start + static_cast<T>(i) * h);
  }
  return h * result / static_cast<T>(3);
}

// Bilinear interpolation on a 2D grid.
// z is stored in row-major order: z[j + ny*i] = f(x[i], y[j]).
template <typename T>
T interpolate2d(const std::vector<T>& x, const std::vector<T>& y, const std::vector<T>& z, T xi, T yj) {
  const size_t nx = x.size();
  const size_t ny = y.size();
  std::vector<double> za(nx * ny);

  gsl_spline2d* spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
  gsl_interp_accel* xacc = gsl_interp_accel_alloc();
  gsl_interp_accel* yacc = gsl_interp_accel_alloc();

  for (size_t i = 0; i < nx; ++i)
    for (size_t j = 0; j < ny; ++j) gsl_spline2d_set(spline, za.data(), i, j, static_cast<double>(z.at(j + ny * i)));

  gsl_spline2d_init(spline, x.data(), y.data(), za.data(), nx, ny);
  const T result = static_cast<T>(gsl_spline2d_eval(spline, xi, yj, xacc, yacc));

  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  return result;
}

}  // namespace Numeric
}  // namespace CRAMS

#endif  // INCLUDE_NUMERIC_H
