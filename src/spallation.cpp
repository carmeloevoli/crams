#include "spallation.h"
#include "utilities.h"

SpallationXsecs::SpallationXsecs(const PID& fragment, const double& Be_xsecs_norm, bool doError) :
		_fragment(fragment), _Be_xsecs_norm(Be_xsecs_norm), _doError(doError) {
	read_table();
#ifdef DEBUG
	std::cout << "read " << _table.size() << " cross-sections for " << _fragment << "\n";
#endif
	double T = _T_min;
	for (size_t i = 0; i < _T_size; ++i) {
		_T.push_back(T);
		T *= _T_ratio;
	}
}

SpallationXsecs::~SpallationXsecs() {
}

double SpallationXsecs::get_error() const {
	const double error = 0.3;
	std::random_device rd { };
	std::mt19937 gen { rd() };
	std::normal_distribution<> d { 1., error };
	double value = d(gen);
	return std::fabs(value);
}

double SpallationXsecs::get_ISM(const PID& projectile, const double& T) const {
	double value = 0;
	auto it = _table.find(projectile);
	if (it != _table.end()) {
		double T_now = std::min(T, _T.back());
		value = LinearInterpolatorLog(_T, it->second, T_now);
		auto it_error = _xsec_error.find(projectile);
		value *= (_doError) ? it_error->second : 1.;
	}
	if (_fragment.get_Z() == 4)
		value *= _Be_xsecs_norm;
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
			x.reserve(_T_size);
			for (size_t i = 0; i < _T_size; ++i) {
				inf >> x_temp;
				x.emplace_back(std::max(x_temp, 1e-10) * cgs::mbarn);
			}
			if (PID(Z_frag, A_frag) == _fragment) {
				_table[PID(Z_proj, A_proj)] = x;
				_xsec_error[PID(Z_proj, A_proj)] = get_error();
			}
		}
		inf.close();
	}
}

