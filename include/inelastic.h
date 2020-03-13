#ifndef INCLUDE_SPALLATION_H_
#define INCLUDE_SPALLATION_H_

#include <cmath>
#include <vector>

#include "cgs.h"
#include "params.h"
#include "pid.h"

double sigma_pp(const double& T);

double sigma_ST(const double& T);

class InelasticXsec {
public:
	InelasticXsec();
	InelasticXsec(const PID& proj, bool doError = false);
	virtual ~InelasticXsec();
	virtual double get_ISM(const double& T) const = 0;

protected:
	PID _proj;
	double _renorm_factor = 1.;

private:
	double get_error() const;
	bool _doError = false;
};

class InXsecTripathi99: public InelasticXsec {
public:
	InXsecTripathi99();
	InXsecTripathi99(const PID& proj, bool doError = false);
	double get_ISM(const double& T) const override;
};

class InXsecCROSEC: public InelasticXsec {
public:

	InXsecCROSEC();
	InXsecCROSEC(const PID& proj, bool doError = false);
	double get_ISM(const double& T) const override;

protected:
	void read_table(const std::string& filename);
	double get_H(const double& T) const;

protected:
	std::vector<double> _table;
	std::vector<double> _T;
	double _T_min = 0.1 * cgs::GeV;
	size_t _T_size = 100;
	double _T_ratio = 1.1;
};

//class CROSEC: public InelasticXsecTable {
//public:
//	CROSEC(const PID& projectile) :
//			InelasticXsecTable(projectile) {
//		read_table("data/xsecs_total_CROSEC_0.1_100_1.1.txt");
//	}
//};
//
//

#endif /* INCLUDE_SPALLATION_H_ */
