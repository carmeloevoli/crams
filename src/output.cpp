#include "output.h"
#include "utilities.h"

OutputManager::OutputManager(const std::vector<Particle>& particles) :
		_particles(particles) {
}

OutputManager::~OutputManager() {
}

void OutputManager::dump_spectra(double R_min, double R_max, size_t R_size) const {
	auto ptr_H1 = find(_particles.begin(), _particles.end(), Particle(H1, 0));
	auto ptr_B10 = find(_particles.begin(), _particles.end(), Particle(B10, 0));
	auto ptr_B11 = find(_particles.begin(), _particles.end(), Particle(B11, 0));
	auto ptr_C12 = find(_particles.begin(), _particles.end(), Particle(C12, 0));
	auto ptr_O16 = find(_particles.begin(), _particles.end(), Particle(O16, 0));
	auto R = LogAxis(R_min, R_max, R_size);
	std::ofstream outfile("spectra.txt");
	outfile << std::scientific;
	double units = 1. / (cgs::GeV * pow2(cgs::meter) * cgs::sec);
	for (size_t i = 0; i < R_size; ++i) {
		double R_ = R.at(i);
		outfile << R_ / cgs::GeV << "\t";
		outfile << ptr_H1->get_I_R_TOA(R_, 0.7 * cgs::GeV) / units << "\t";
		outfile << ptr_B10->get_I_R_TOA(R_, 0.7 * cgs::GeV) / units << "\t";
		outfile << ptr_B11->get_I_R_TOA(R_, 0.7 * cgs::GeV) / units << "\t";
		outfile << ptr_C12->get_I_R_TOA(R_, 0.7 * cgs::GeV) / units << "\t";
		outfile << ptr_O16->get_I_R_TOA(R_, 0.7 * cgs::GeV) / units << "\t";
		outfile << "\n";
	}
	outfile.close();
}

void OutputManager::dump_ratio(double R_min, double R_max, size_t R_size) const {
	auto ptr_B10 = find(_particles.begin(), _particles.end(), Particle(B10, 0));
	auto ptr_B11 = find(_particles.begin(), _particles.end(), Particle(B11, 0));
	auto ptr_C12 = find(_particles.begin(), _particles.end(), Particle(C12, 0));
	auto ptr_O16 = find(_particles.begin(), _particles.end(), Particle(O16, 0));
	auto R = LogAxis(R_min, R_max, R_size);
	std::ofstream outfile("ratios.txt");
	outfile << std::scientific;
	for (size_t i = 0; i < R_size; ++i) {
		double R_ = R.at(i);
		outfile << R_ / cgs::GeV << "\t";
		double B = ptr_B10->get_I_R_TOA(R_, 0.7 * cgs::GeV)
				+ ptr_B11->get_I_R_TOA(R_, 0.7 * cgs::GeV);
		double C = ptr_C12->get_I_R_TOA(R_, 0.7 * cgs::GeV);
		double O = ptr_O16->get_I_R_TOA(R_, 0.7 * cgs::GeV);
		outfile << B / C << "\t";
		outfile << C / O << "\t";
		outfile << "\n";
	}
	outfile.close();
}
