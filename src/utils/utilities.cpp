#include "crams/utils/utilities.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>

#include "crams/core/cgs.h"

namespace CRAMS {
namespace Utilities {

namespace {
constexpr size_t GSL_LIMIT = 1000;
constexpr double GSL_EPSREL = 1e-6;
}  // namespace

double T2beta(double T) {
  if (T < 0.0) throw std::invalid_argument("T must be positive");
  const double beta = std::sqrt(T * (T + 2. * CGS::protonMassC2)) / (T + CGS::protonMassC2);
  return beta;
}

double T2gamma(double T) {
  if (T < 0.0) throw std::invalid_argument("T must be positive");
  return (T + CGS::protonMassC2) / CGS::protonMassC2;
}

double T2pc(double T, const PID& pid) {
  if (T < 0.0) throw std::invalid_argument("T must be positive");
  return std::sqrt(T * (T + 2. * CGS::protonMassC2)) * (double)pid.getA();
}

double R2T(double R, const PID& pid) {
  if (R < 0.0) throw std::invalid_argument("R must be positive");
  const double mpSquared = pow2(CGS::protonMassC2);
  const double ZOverASquared = pow2(pid.getZoverA());
  return std::sqrt(pow2(R) * ZOverASquared + mpSquared) - CGS::protonMassC2;
}

double computeRandomFactor(double variance) {
  thread_local static std::mt19937 gen{std::random_device{}()};
  std::normal_distribution<> d{1., variance};
  return std::fabs(d(gen));
}

std::vector<double> LinAxis(double min, double max, size_t size) {
  if (!(min < max)) throw std::invalid_argument("min must be smaller than max");
  if (!(size > 1)) throw std::invalid_argument("size must be larger than 1");

  const double dx = (max - min) / (double)(size - 1);
  std::vector<double> v(size);
  for (size_t i = 0; i < size; ++i) v[i] = min + dx * i;
  return v;
}

std::vector<double> LogAxis(double min, double max, size_t size) {
  if (!(min < max)) throw std::invalid_argument("min must be smaller than max");
  if (!(size > 1)) throw std::invalid_argument("size must be larger than 1");

  const double log_min = std::log(min);
  const double log_step = std::log(max / min) / (double)(size - 1);
  std::vector<double> v(size);
  for (size_t i = 0; i < size; ++i) v[i] = std::exp(log_min + i * log_step);
  return v;
}

bool isGoodAndPositive(const std::vector<double>& v) {
  return std::none_of(v.begin(), v.end(), [](double d) {
    return std::isnan(d) || std::isinf(d) || d < 0.;
  });
}

bool fileExists(const std::string& filename) {
  std::ifstream f(filename);
  return f.good();
}

std::string simplifyKey(const std::string& key) {
  std::string value;
  for (char c : key) {
    if (c != '_' && c != ' ') value += std::tolower(c);
  }
  return value;
}

static std::vector<std::string> splitrow(const std::string& s) {
  std::vector<std::string> result;
  std::istringstream iss(s);
  std::string token;
  while (iss >> token) result.push_back(token);
  return result;
}

std::vector<double> loadColumn(const std::string& filename, size_t useCol, size_t nHeaderLines) {
  std::ifstream file(filename);
  if (!file.is_open()) throw std::runtime_error("Cannot open file: " + filename);

  std::vector<double> v;
  std::string line;
  size_t lineNum = 0;
  while (std::getline(file, line)) {
    if (lineNum++ < nHeaderLines) continue;
    const auto items = splitrow(line);
    if (!items.empty() && useCol < items.size()) {
      v.push_back(std::stod(items[useCol]));
    }
  }
  return v;
}

bool inRange(double x, std::pair<double, double> range) {
  return x >= range.first && x <= range.second;
}

}  // namespace Utilities
}  // namespace CRAMS
