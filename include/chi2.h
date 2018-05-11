#ifndef INCLUDE_CHI2_H_
#define INCLUDE_CHI2_H_

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include "mks.h"
#include "pid.h"

struct data_point {
	double rigidity;
	double flux;
	double flux_err_low;
	double flux_err_high;
};

class Chi2 {
public:
	Chi2();
	Chi2(const std::string& filename);
	~Chi2();
	void read_datafile(const std::string& filename);
	void calculate_chi2(const PID& pid, const std::vector<double>& T, const std::vector<double>& flux, const double& modulation_potential);

	inline double getChi2() const {
		return chi2;
	}

private:
	double chi2 = 0;
	std::vector<data_point> data;
};

#endif /* INCLUDE_CHI2_H_ */
