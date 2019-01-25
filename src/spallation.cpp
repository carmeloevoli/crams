#include "spallation.h"
#include "utilities.h"

SpallationXsecs::SpallationXsecs(const PID& fragment, bool doError) :
		_fragment(fragment), _doError(doError) {
	read_table();
	std::cout << "read " << _table.size() << " cross-sections for " << _fragment << "\n";
	double T = _T_min;
	for (size_t i = 0; i < _T_size; ++i) {
		_T.push_back(T);
		T *= _T_ratio;
	}
}

SpallationXsecs::~SpallationXsecs() {
}

double SpallationXsecs::get_error(const PID& projectile) const {
	double value = 1;
	auto it = _xsec_error.find(projectile);
	if (it != _xsec_error.end()) {
		double error = it->second;
		std::random_device rd { };
		std::mt19937 gen { rd() };
		std::normal_distribution<> d { 1., error };
		value = d(gen);
	}
	return std::fabs(value);
}

double SpallationXsecs::get_ISM(const PID& projectile, const double& T) const {
	double value = 0;
	auto it = _table.find(projectile);
	if (it != _table.end()) {
		double T_now = std::min(T, _T.back());
		value = LinearInterpolator(_T, it->second, T_now);
	}
	double error = (_doError) ? get_error(projectile) : 1.;
	return error * value * (1. + cgs::K_He * cgs::f_He) / (1. + cgs::f_He);
}

void SpallationXsecs::read_table() {
	std::ifstream inf(_table_filename.c_str());
	if (!inf.is_open()) {
		std::cout << "filename " << _table_filename << " cannot be open.\n";
		exit(2);
	} else {
		int Z_frag, A_frag, Z_proj, A_proj;
		double x_temp;
		while (inf) {
			inf >> Z_frag >> A_frag >> Z_proj >> A_proj;
			std::vector<double> x;
			x.reserve(_T_size);
			for (size_t i = 0; i < _T_size; ++i) {
				inf >> x_temp;
				x.emplace_back(x_temp * cgs::mbarn);
			}
			if (PID(Z_frag, A_frag) == _fragment) {
				_table[PID(Z_proj, A_proj)] = x;
				_xsec_error[PID(Z_proj, A_proj)] = 0.2;
			}
		}
		inf.close();
	}
}

