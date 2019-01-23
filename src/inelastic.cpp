#include "inelastic.h"
#include "utilities.h"

InelasticXsec::InelasticXsec() {
}

InelasticXsec::InelasticXsec(const PID& pid) {
	A = pid.get_A();
	Z = pid.get_Z();
	table = Tripathi99(pid);
}

InelasticXsec::~InelasticXsec() {
	//std::cout << "delete inelastic sigma for particle " << A << " " << Z << "\n";
}

double InelasticXsec::get_ISM(const double& T) const {
	//double sigma = (Z == 1) ? sigma_pp(T) : sigma_ST(T);
	double sigma = (Z == 1) ? sigma_pp(T) : table.get_H(T);
	sigma *= (1. + cgs::K_He * cgs::f_He) / (1. + cgs::f_He);
	return std::max(sigma, 1e-10 * cgs::mbarn);
}

double sigma_pp(const double& T) {
	constexpr double E_threshold = 0.2797 * cgs::GeV;
	double value = 0;
	double x = T / E_threshold;
	if (x > 1) {
		value = 30.7 - 0.96 * log(x) + 0.18 * pow2(log(x));
		value *= pow3(1 - pow(x, -1.9));
	}
	return value * cgs::mbarn;
}

double InelasticXsec::sigma_ST(const double& T) const {
	double value = 45. * std::pow((double) A, 0.7);
	value *= 1. + 0.016 * std::sin(5.3 - 2.63 * std::log(A));
	double T_MeV = T / cgs::MeV;
	value *= 1. - 0.62 * std::exp(-T_MeV / 200.) * std::sin(10.9 * pow(T_MeV, -0.28));
	return value * cgs::mbarn;
}

InelasticXsecTable::InelasticXsecTable() {
}

InelasticXsecTable::InelasticXsecTable(const PID& projectile) :
		_projectile(projectile) {
	double T = _T_min;
	for (size_t i = 0; i < _T_size; ++i) {
		_T.push_back(T);
		T *= _T_ratio;
	}
}

InelasticXsecTable::~InelasticXsecTable() {
}

double InelasticXsecTable::get_H(const double& T) const {
	double T_now = std::min(T, _T.back());
	return LinearInterpolator(_T, _table, T_now);
}

void InelasticXsecTable::read_table(const std::string& filename) {
	std::ifstream inf(filename.c_str());
	if (!inf.is_open()) {
		std::cout << "filename " << filename << " cannot be open.\n";
		exit(2);
	} else {
		int Z_proj, A_proj;
		double x_temp;
		while (inf) {
			inf >> Z_proj >> A_proj;
			std::vector<double> x;
			x.reserve(_T_size);
			for (size_t i = 0; i < _T_size; ++i) {
				inf >> x_temp;
				x.emplace_back(x_temp * cgs::mbarn);
			}
			if (PID(Z_proj, A_proj) == _projectile)
				copy(x.begin(), x.end(), back_inserter(_table));
		}
		inf.close();
	}
	std::cout << "inelastic table read with " << _table.size() << " points.\n";
}
