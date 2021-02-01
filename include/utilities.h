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
double LinearInterpolator(const std::vector<double>& x, const std::vector<double>& y, const double& x_new);
double LinearInterpolatorLog(const std::vector<double>& x, const std::vector<double>& y, const double& x_new);
bool isGoodAndPositive(const std::vector<double>& v);
bool fileExists(const std::string& filename);
std::string simplifyKey(const std::string& key);
std::vector<std::string> split(const std::string& str, const std::string& delim);

}  // namespace Utilities
}  // namespace CRAMS

#endif /* INCLUDE_UTILITIES_H_ */
