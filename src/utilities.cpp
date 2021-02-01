#include "utilities.h"

#include <gsl/gsl_integration.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

#include "cgs.h"

#define LIMIT 1000

namespace CRAMS {
namespace Utilities {

double T2beta(const double& T) {
  if (T < 0.0) throw std::invalid_argument("T must be positive");

  double beta = std::sqrt(T * (T + 2. * CGS::protonMassC2));
  beta /= T + CGS::protonMassC2;
  return beta;
}

double T2gamma(const double& T) {
  if (T < 0.0) throw std::invalid_argument("T must be positive");

  double gamma = T + CGS::protonMassC2;
  gamma /= CGS::protonMassC2;
  return gamma;
}

double T2pc(const double& T, const PID& pid) {
  if (T < 0.0) throw std::invalid_argument("T must be positive");

  double pc = std::sqrt(T * (T + 2. * CGS::protonMassC2));
  pc *= (double)pid.getA();
  return pc;
}

double R2T(const double& R, const PID& pid) {
  if (R < 0.0) throw std::invalid_argument("R must be positive");

  constexpr double mpSquared = pow2(CGS::protonMassC2);
  const double ZOverASquared = pow2(pid.getZoverA());
  const double ESquared = pow2(R) * ZOverASquared + mpSquared;
  const double T = std::sqrt(ESquared) - CGS::protonMassC2;
  return T;
}

double computeRandomFactor(const double& variance) {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{1., variance};
  double value = d(gen);
  return std::fabs(value);
}

double gslGammaIntegrand(double x, void* params) {
  const double alpha = *(double*)params;
  double f = std::pow(x, 2. - alpha);
  f *= std::sqrt(x * x + 1) - 1;
  return f;
}

double gslXGammaIntegrand(double y, void* params) {
  const double alpha = *(double*)params;
  double x = std::exp(y);
  double f = std::pow(x, 2. - alpha);
  f *= std::sqrt(x * x + 1) - 1;
  return x * f;
}

double GammaIntegral(double slope) {
  if (!(slope > 4.0 && slope < 5.0)) throw std::invalid_argument("slope must be greater than 4 and smaller than 5");

  gsl_integration_workspace* w = gsl_integration_workspace_alloc(LIMIT);
  double result, error;
  gsl_function F;
  F.params = &slope;

  const double EPSREL = 1e-6;

  if (slope < 4.1) {
    F.function = &gslGammaIntegrand;
    gsl_integration_qagiu(&F, 0, 0, EPSREL, LIMIT, w, &result, &error);
  } else {
    F.function = &gslXGammaIntegrand;
    gsl_integration_qag(&F, std::log(1e-5), std::log(1e10), 0, EPSREL, LIMIT, 4, w, &result, &error);
  }
  gsl_integration_workspace_free(w);
  return 4. * M_PI * result;
}

std::vector<double> LinAxis(const double& min, const double& max, const size_t& size) {
  if (!(min < max)) throw std::invalid_argument("min must be smaller than max");
  if (!(size > 1)) throw std::invalid_argument("size must be larger than 1");

  const double dx = (max - min) / (double)(size - 1);
  std::vector<double> v(size);
  for (size_t i = 0; i < size; ++i) {
    const auto value = min + dx * i;
    v[i] = value;
  }
  return v;
}

std::vector<double> LogAxis(const double& min, const double& max, const size_t& size) {
  if (!(min < max)) throw std::invalid_argument("min must be smaller than max");
  if (!(size > 1)) throw std::invalid_argument("size must be larger than 1");

  const double delta_log = std::exp(std::log(max / min) / (size - 1));
  std::vector<double> v(size);
  for (size_t i = 0; i < size; ++i) {
    const auto value = std::exp(std::log(min) + (double)i * std::log(delta_log));
    v[i] = value;
  }
  return v;
}

double LinearInterpolator(const std::vector<double>& x, const std::vector<double>& y, const double& x_new) {
  if (x_new < x.front() || x_new > x.back())
    throw std::invalid_argument("x_new out of vector range in LinearInterpolator");

  size_t const i = std::lower_bound(x.begin(), x.end(), x_new) - x.begin();
  const auto t = (x_new - x.at(i - 1)) / (x.at(i) - x.at(i - 1));
  return y.at(i - 1) * (1. - t) + y.at(i) * t;
}

// double LinearInterpolatorLog(const std::vector<double>& x, const std::vector<double>& y,
//		const double& x_new) {
//	auto value = double();
//	if (x_new >= x.front() && x_new <= x.back()) {
//		size_t const i = std::lower_bound(x.begin(), x.end(), x_new) - x.begin();
//		double t = std::log(x_new) - std::log(x.at(i - 1));
//		t /= std::log(x.at(i)) - std::log(x.at(i - 1));
//		value = std::log(y.at(i - 1)) * (1. - t) + std::log(y.at(i)) * t;
//		value = std::exp(value);
//	}
//	return value;
//}

size_t getLowerIndex(const std::vector<double>& v, double x) {  // TODO speed up this function?
  assert(x >= v.front());
  size_t i = 0;
  while (v.at(i + 1) < x) i++;
  return i;
}

double LinearInterpolatorLog(const std::vector<double>& x, const std::vector<double>& y, const double& x_new) {
  if (x_new < x.front() || x_new > x.back())
    throw std::invalid_argument("x_new out of vector range in LinearInterpolatorLog");

  size_t const i = getLowerIndex(x, x_new);  // TODO uniformize with LinearInterpolator
  double t = std::log(x_new) - std::log(x.at(i));
  t /= std::log(x.at(i + 1)) - std::log(x.at(i));
  double v = std::log(y.at(i)) * (1. - t) + std::log(y.at(i + 1)) * t;
  return std::exp(v);
}

bool isGoodAndPositive(const std::vector<double>& v) {
  bool isBad = std::any_of(v.begin(), v.end(), [](double d) { return std::isnan(d); });
  bool isNegative = std::any_of(v.begin(), v.end(), [](double d) { return (d < 0.); });
  return !isBad && !isNegative;
}

bool fileExists(const std::string& filename) {
  std::ifstream f(filename.c_str());
  return f.good();
}

std::string simplifyKey(const std::string& key) {
  std::string value;
  for (auto it = key.begin(); it != key.end(); ++it) {
    if (*it != '_' && *it != ' ') {
      value += tolower(*it);
    }
  }
  return value;
}

std::vector<std::string> split(const std::string& str, const std::string& delim) {
  if (str == "" || delim == "") throw std::invalid_argument("str or delim cannot be empty");

  std::vector<std::string> tokens;
  size_t prev = 0, pos = 0;
  do {
    pos = str.find(delim, prev);
    if (pos == std::string::npos) pos = str.length();
    auto token = str.substr(prev, pos - prev);
    if (!token.empty()) tokens.push_back(token);
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());
  return tokens;
}

#undef LIMIT

}  // namespace Utilities
}  // namespace CRAMS