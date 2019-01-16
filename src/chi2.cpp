#include "chi2.h"

#define max_num_of_char_in_a_line 512

Chi2::Chi2() {
}

Chi2::Chi2(const std::string& filename) {
	read_datafile(filename);
}

Chi2::~Chi2() {
}

void Chi2::read_datafile(const std::string& filename) {
	std::cout << "reading data from " << filename << "... ";
	std::ifstream file_to_read(filename.c_str());
	if (file_to_read.is_open()) {
		for (int i = 0; i < 2; i++)
		file_to_read.ignore(max_num_of_char_in_a_line, '\n');
		double read_values[4];
		while (!file_to_read.eof()) {
			file_to_read >> read_values[0] >> read_values[1] >> read_values[2] >> read_values[3];
			data_point point;
//			point.rigidity = read_values[0] * GeV;
//			point.flux = read_values[1] * (1. / GeV / m2 / sec / sr);
//			point.flux_err_low = read_values[2] * (1. / GeV / m2 / sec / sr);
//			point.flux_err_high = read_values[3] * (1. / GeV / m2 / sec / sr);
//			data.push_back(point);
		}
	}
	std::cout << " with size : " << data.size() << "\n";
	file_to_read.close();
}

/*void Chi2::read_model(string fitsname, double modpotential)
 {
 read_spectrum(fitsname,1,1,&protons.E,&protons.dNdE);
 modulated_spectrum(&protons.E,&protons.dNdE,&protons.modulated,1,1,modpotential);

 return;
 }*/

void Chi2::calculate_chi2(const PID& pid, const std::vector<double>& T, const std::vector<double>& flux, const double& modulation_potential) {
	auto T_size = T.size();
	double log_T[T_size];
	double log_I[T_size];

	for (int i = 0; i < T_size; ++i) {
		log_T[i] = std::log(T.at(i));
		log_I[i] = (flux.at(i) > 0) ? std::log(flux.at(i)) : -100.;
	}

	double Phi = pid.get_Z_over_A() * modulation_potential;

	gsl_interp_accel * acc = gsl_interp_accel_alloc();
	gsl_spline * spline = gsl_spline_alloc(gsl_interp_cspline, T_size);
	gsl_spline_init(spline, log_T, log_I, T_size);

	chi2 = 0.0;

	for (auto it = data.begin(); it != data.end(); it++) {
		double R = it->rigidity;

//		std::cout << "R : " << R / GeV << "\t";
//
//		double T_ = std::sqrt(pow2(pid.get_Z_over_A() * R) + pow2(proton_mass_c2)) - proton_mass_c2;
//		const double T_phi = std::min(T_ + Phi, T.back());
//		double factor = pid.get_Z_over_A();
//		factor *= std::sqrt(T_ * (T_ + 2. * proton_mass_c2)) / (T_ + proton_mass_c2);
//		factor *= T_ * (T_ + 2.0 * proton_mass_c2) / T_phi / (T_phi + 2.0 * proton_mass_c2);
//		double I_R_TOA = factor * std::exp(gsl_spline_eval(spline, std::log(T_phi), acc));
//
//		double delta_chi2 = pow2(I_R_TOA - it->flux);
//		delta_chi2 /= (I_R_TOA < it->flux) ? pow2(it->flux_err_low) : pow2(it->flux_err_high);

//		chi2 += delta_chi2;
	}

	chi2 /= (double)data.size();

	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}
