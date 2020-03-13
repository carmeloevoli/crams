#ifndef _INCLUDE_SPALLATION_H_
#define _INCLUDE_SPALLATION_H_

#include <map>
#include <random>
#include <string>
#include <vector>

#include "cgs.h"
#include "pid.h"

#define GALPROP2017

class SpallationXsecs {
public:
	SpallationXsecs(const PID& fragment, const double& Be_xsecs_norm, bool doError = false);
	virtual ~SpallationXsecs();

	double get_ISM(const PID& projectile, const double& T) const;

protected:
	void read_table();
	double get_error() const;

protected:
	bool _doError = false;
	double _renorm_factor = 1;
	PID _fragment;
	std::map<PID, double> _xsec_error;
	std::map<PID, std::vector<double> > _table;
	std::vector<double> _T;
	double _Be_xsecs_norm = 1;

#ifdef EVOLI2018
	std::string _table_filename = "data/xsecs_0.1_100_1.1.txt";
	double _T_min = 0.1 * cgs::GeV;
	size_t _T_size = 100;
	double _T_ratio = 1.1;
#endif

#ifdef GALPROP2017
	std::string _table_filename = "data/sigProdGALPROP17_OPT12.txt";
	double _T_min = 0.01 * cgs::GeV;
	double _T_max = 10. * cgs::GeV;
	size_t _T_size = 41;
	double _T_ratio = std::pow(_T_max / _T_min, 1. / (_T_size - 1));
#endif

#ifdef USINE
	std::string _table_filename = "data/sigProdWebber03+Coste12.txt"
	double _T_min = 0.01 * cgs::GeV;
	double _T_max = 10. * cgs::GeV;
	size_t _T_size = 41;
	double _T_ratio = std::pow(_T_max/_T_min, 1. / (_T_size - 1));
#endif

};

#endif /* INCLUDE_SPALLATION_H_ */
