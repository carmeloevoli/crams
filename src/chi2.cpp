#include "chi2.h"

#define max_num_of_char_in_a_line 512

Chi2::Chi2() {
}

Chi2::Chi2(const Particles& particles, const double& phi, const std::string& filename) :
		_particles(particles) {
	read_datafile(filename);
	set_phi(phi);
}

Chi2::~Chi2() {
}

void Chi2::read_datafile(const std::string& filename) {
	std::cout << "reading data from " << filename << "... ";
	std::ifstream file_to_read(filename.c_str());
	if (file_to_read.is_open()) {
		file_to_read.ignore(max_num_of_char_in_a_line, '\n');
		double values[4];
		constexpr double units = 1. / (cgs::GeV * cgs::m2 * cgs::sec * cgs::sr);
		while (!file_to_read.eof()) {
			file_to_read >> values[0] >> values[1] >> values[2] >> values[3];
			data_point point;
			point.R = values[0] * cgs::GeV;
			point.F = values[1] * units;
			point.F_err_low = values[2] * units;
			point.F_err_high = values[3] * units;
			_data.push_back(point);
		}
	}
	std::cout << " with size : " << _data.size() << "\n";
	file_to_read.close();
}

double Chi2::compute_chi2(const double& R_min) {
	double chi2 = 0.0;
	size_t ndata = 0;
	for (auto it = _data.begin(); it != _data.end(); it++) {
		if (it->R > R_min) {
			double I_R_TOA = get_I_R_TOA(it->R, _phi);
			double delta_chi2 = pow2(I_R_TOA - it->F);
			delta_chi2 /= (I_R_TOA < it->R) ? pow2(it->F_err_low) : pow2(it->F_err_high);
			chi2 += delta_chi2;
			ndata++;
		}
	}
	return chi2 / (double) ndata;
}

double Chi2_C::get_I_R_TOA(const double& R, const double& phi) {
	double value = (ptr_C12.isPresent) ? ptr_C12.it->I_R_TOA(R, phi) : 0.;
	value += (ptr_C13.isPresent) ? ptr_C13.it->I_R_TOA(R, phi) : 0.;
	value += (ptr_C14.isPresent) ? ptr_C14.it->I_R_TOA(R, phi) : 0.;
	return value;
}

double Chi2_O::get_I_R_TOA(const double& R, const double& phi) {
	double value = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, phi) : 0.;
	value += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, phi) : 0.;
	value += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, phi) : 0.;
	return value;
}
