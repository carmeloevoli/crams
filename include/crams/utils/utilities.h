#ifndef CRAMS_UTILS_UTILITIES_H_
#define CRAMS_UTILS_UTILITIES_H_

#include <string>
#include <vector>

#include "crams/core/pid.h"

namespace CRAMS {
namespace Utilities {

template <typename T>
constexpr T pow2(T x) { return x * x; }

template <typename T>
constexpr T pow3(T x) { return x * x * x; }

template <typename T>
constexpr T pow4(T x) { return x * x * x * x; }

double T2gamma(double T);
double T2beta(double T);
double T2pc(double T, const PID& pid);
double R2T(double R, const PID& pid);

double computeRandomFactor(double variance);
std::vector<double> LinAxis(double min, double max, size_t size);
std::vector<double> LogAxis(double min, double max, size_t size);
bool isGoodAndPositive(const std::vector<double>& v);
bool fileExists(const std::string& filename);
std::string simplifyKey(const std::string& key);
std::vector<double> loadColumn(const std::string& filename, size_t useCol, size_t nHeaderLines);
bool inRange(double x, std::pair<double, double> range);

}  // namespace Utilities
}  // namespace CRAMS

#endif  // CRAMS_UTILS_UTILITIES_H_
