#include "spallation.h"
#include "utilities.h"

SpallationXsecs::SpallationXsecs(const PID& fragment) :
		_fragment(fragment) {
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

double SpallationXsecs::get(const PID& projectile, const double& T) const {
	double value = 0;
	auto it = _table.find(projectile);
	if (it != _table.end()) {
		double T_now = std::min(T, _T.back());
		return LinearInterpolator(_T, it->second, T_now);
	}
	return value * (1. + cgs::K_He * cgs::f_He) / (1. + cgs::f_He);
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
			x.reserve(100);
			for (size_t i = 0; i < 100; ++i) {
				inf >> x_temp;
				x.emplace_back(x_temp * cgs::mbarn);
			}
			if (PID(Z_frag, A_frag) == _fragment) {
				_table[PID(Z_proj, A_proj)] = x;
			}
		}
		inf.close();
	}
}

