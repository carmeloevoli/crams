#ifndef INCLUDE_GSL_H
#define INCLUDE_GSL_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include <functional>

namespace GSL {

template <typename T>
size_t getLowerIndex(const std::vector<T> &v, T x) {  // TODO speed up this function?
  assert(x >= v.front());
  size_t i = 0;
  while (v.at(i + 1) < x) i++;
  return i;
}

template <typename T>
T LinearInterpolator(const std::vector<T> &x, const std::vector<T> &y, T x_new) {
  if (x_new < x.front() || x_new > x.back())
    throw std::invalid_argument("x_new out of vector range in LinearInterpolator");

  size_t const i = getLowerIndex(x, x_new);  // std::lower_bound(x.begin(), x.end(), x_new) - x.begin();
  const auto t = (x_new - x[i - 1]) / (x[i] - x[i - 1]);
  return y[i - 1] * (1. - t) + y[i] * t;
}

template <typename T>
T LinearInterpolatorLog(const std::vector<T> &x, const std::vector<T> &y, T x_new) {
  if (x_new < x.front() || x_new > x.back())
    throw std::invalid_argument("x_new out of vector range in LinearInterpolatorLog");

  size_t const i = getLowerIndex(x, x_new);  // std::lower_bound(x.begin(), x.end(), x_new) - x.begin();
  double t = std::log(x_new) - std::log(x[i]);
  t /= std::log(x[i + 1]) - std::log(x[i]);
  double v = std::log(y[i]) * (1. - t) + std::log(y[i + 1]) * t;
  return std::exp(v);
}

template <typename T>
T QAGIntegration(std::function<T(T)> f, T start, T stop, int LIMIT, double rel_error = 1e-4) {
  double a = static_cast<double>(start);
  double b = static_cast<double>(stop);
  double abs_error = 0.0;  // disabled
  int key = GSL_INTEG_GAUSS31;
  double result;
  double error;

  gsl_function F;
  F.function = [](double x, void *vf) -> double {
    auto &func = *static_cast<std::function<double(double)> *>(vf);
    return func(x);
  };
  F.params = &f;

  gsl_integration_workspace *workspace_ptr = gsl_integration_workspace_alloc(LIMIT);
  gsl_integration_qag(&F, a, b, abs_error, rel_error, LIMIT, key, workspace_ptr, &result, &error);
  gsl_integration_workspace_free(workspace_ptr);

  return T(result);
}

template <typename T>
T QAGSIntegration(std::function<T(T)> f, T start, T stop, int LIMIT, double rel_error = 1e-4) {
  double a = static_cast<double>(start);
  double b = static_cast<double>(stop);
  double abs_error = 0.0;  // disabled
  double result;
  double error;

  gsl_function F;
  F.function = [](double x, void *vf) -> double {
    auto &func = *static_cast<std::function<double(double)> *>(vf);
    return func(x);
  };
  F.params = &f;

  gsl_integration_workspace *workspace_ptr = gsl_integration_workspace_alloc(LIMIT);
  gsl_integration_qags(&F, a, b, abs_error, rel_error, LIMIT, workspace_ptr, &result, &error);
  gsl_integration_workspace_free(workspace_ptr);

  return T(result);
}

template <typename T>
T simpsonIntegration(std::function<T(T)> f, T start, T stop, int N = 100) {
  const T a = start;
  const T b = stop;

  const T h = (b - a) / N;
  const T XI0 = f(a) + f(b);

  T XI1 = 0, XI2 = 0;

  for (int i = 1; i < N; ++i) {
    const T X = a + i * h;
    if (i % 2 == 0)
      XI2 = XI2 + f(X);
    else
      XI1 = XI1 + f(X);
  }

  return h * (XI0 + 2 * XI2 + 4 * XI1) / 3.0;
}

}  // namespace GSL

#endif