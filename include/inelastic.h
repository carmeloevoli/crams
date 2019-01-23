#ifndef INCLUDE_SPALLATION_H_
#define INCLUDE_SPALLATION_H_

#include <cmath>
#include <vector>

#include "cgs.h"
#include "params.h"
#include "pid.h"

double sigma_pp(const double& T);

class InelasticXsecTable {
public:
	InelasticXsecTable();
	InelasticXsecTable(const PID& projectile);
	virtual ~InelasticXsecTable();

	double get_H(const double& T) const;

protected:
	void read_table(const std::string& filename);

protected:
	PID _projectile;
	std::vector<double> _table;
	std::vector<double> _T;
	double _T_min = 0.1 * cgs::GeV;
	size_t _T_size = 100;
	double _T_ratio = 1.1;
};

class Tripathi99: public InelasticXsecTable {
public:
	Tripathi99(const PID& projectile) :
			InelasticXsecTable(projectile) {
		read_table("data/xsecs_total_Tripathi99_0.1_100_1.1.txt");
	}
};

class CROSEC: public InelasticXsecTable {
public:
	CROSEC(const PID& projectile) :
			InelasticXsecTable(projectile) {
		read_table("data/xsecs_total_CROSEC_0.1_100_1.1.txt");
	}
};

class InelasticXsec {
public:
	InelasticXsec();
	InelasticXsec(const PID& pid);
	virtual ~InelasticXsec();
	double get_ISM(const double& T) const;

protected:
	int A = 0;
	int Z = 0;

private:
	double sigma_ST(const double& T) const;
	InelasticXsecTable table;
};

#endif /* INCLUDE_SPALLATION_H_ */
