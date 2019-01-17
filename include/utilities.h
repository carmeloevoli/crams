#ifndef INCLUDE_UTILITIES_H_
#define INCLUDE_UTILITIES_H_

#include <vector>

#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define pow4(A) ((A)*(A)*(A)*(A))

double gamma_func(const double& T);
double beta_func(const double& T);
double pc_func(const int& A, const double& T);
double Gamma_Integral(double slope);
std::vector<double> LinAxis(const double& min, const double& max, const size_t& size);
std::vector<double> LogAxis(const double& min, const double& max, const size_t& size);
double LinearInterpolator(const std::vector<double>& x, const std::vector<double>& y, const double& x_new);
double LinearInterpolatorLog(const std::vector<double>& x, const std::vector<double>& y, const double& x_new);

#endif /* INCLUDE_UTILITIES_H_ */
