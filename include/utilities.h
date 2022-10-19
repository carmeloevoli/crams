#ifndef INCLUDE_UTILITIES_H_
#define INCLUDE_UTILITIES_H_

#include <string>
#include <vector>

#include "pid.h"

namespace CRAMS {
namespace Utilities {

#define pow2(A) ((A) * (A))
#define pow3(A) ((A) * (A) * (A))
#define pow4(A) ((A) * (A) * (A) * (A))

double T2gamma(const double& T);
double T2beta(const double& T);
double T2pc(const double& T, const PID& pid);

double computeRandomFactor(const double& variance);
double GammaIntegral(double slope);
std::vector<double> LinAxis(const double& min, const double& max, const size_t& size);
std::vector<double> LogAxis(const double& min, const double& max, const size_t& size);
bool isGoodAndPositive(const std::vector<double>& v);
bool fileExists(const std::string& filename);
std::string simplifyKey(const std::string& key);
std::vector<double> loadColumn(const std::string& filename, size_t useCol, size_t nHeaderLines);
bool inRange(double x, std::pair<double, double> range);

}  // namespace Utilities
}  // namespace CRAMS

#endif /* INCLUDE_UTILITIES_H_ */
