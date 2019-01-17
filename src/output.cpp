#include "output.h"
#include "utilities.h"

OutputManager::OutputManager(const std::vector<Particle>& particles, const double phi) :
		_particles(particles), _phi(phi) {
}

OutputManager::~OutputManager() {
}

void OutputManager::dump_spectra(double R_min, double R_max, size_t R_size) const {
	auto ptr_H1 = find(_particles.begin(), _particles.end(), Particle(H1, 0));
	auto ptr_B10 = find(_particles.begin(), _particles.end(), Particle(B10, 0));
	auto ptr_B11 = find(_particles.begin(), _particles.end(), Particle(B11, 0));
	auto ptr_C12 = find(_particles.begin(), _particles.end(), Particle(C12, 0));
	auto ptr_O16 = find(_particles.begin(), _particles.end(), Particle(O16, 0));
	auto _R = LogAxis(R_min, R_max, R_size);
	std::ofstream outfile("spectra.txt");
	outfile << std::scientific;
	double units = 1. / (cgs::GeV * pow2(cgs::meter) * cgs::sec);
	for (auto& R : _R) {
		outfile << R / cgs::GeV << "\t";
		outfile << ptr_H1->get_I_R_TOA(R, _phi) / units << "\t";
		outfile << ptr_B10->get_I_R_TOA(R, _phi) / units << "\t";
		outfile << ptr_B11->get_I_R_TOA(R, _phi) / units << "\t";
		outfile << ptr_C12->get_I_R_TOA(R, _phi) / units << "\t";
		outfile << ptr_O16->get_I_R_TOA(R, _phi) / units << "\t";
		outfile << "\n";
	}
	outfile.close();
}

void OutputManager::dump_heavy_spectra(double R_min, double R_max, size_t R_size) const {
	auto ptr_Ne20 = find(_particles.begin(), _particles.end(), Particle(Ne20, 0));
	auto ptr_Mg24 = find(_particles.begin(), _particles.end(), Particle(Mg24, 0));
	auto ptr_Si28 = find(_particles.begin(), _particles.end(), Particle(Si28, 0));
	auto ptr_Fe56 = find(_particles.begin(), _particles.end(), Particle(Fe56, 0));
	auto _R = LogAxis(R_min, R_max, R_size);
	std::ofstream outfile("heavy_spectra.txt");
	outfile << std::scientific;
	double units = 1. / (cgs::GeV * pow2(cgs::meter) * cgs::sec);
	for (auto& R : _R) {
		outfile << R / cgs::GeV << "\t";
		outfile << ptr_Ne20->get_I_R_TOA(R, _phi) / units << "\t";
		outfile << ptr_Mg24->get_I_R_TOA(R, _phi) / units << "\t";
		outfile << ptr_Si28->get_I_R_TOA(R, _phi) / units << "\t";
		outfile << ptr_Fe56->get_I_R_TOA(R, _phi) / units << "\t";
		outfile << "\n";
	}
	outfile.close();
}

void OutputManager::dump_ratio(double R_min, double R_max, size_t R_size) const {
	auto ptr_B10 = find(_particles.begin(), _particles.end(), Particle(B10, 0));
	auto ptr_B11 = find(_particles.begin(), _particles.end(), Particle(B11, 0));
	auto ptr_C12 = find(_particles.begin(), _particles.end(), Particle(C12, 0));
	auto ptr_O16 = find(_particles.begin(), _particles.end(), Particle(O16, 0));
	auto _R = LogAxis(R_min, R_max, R_size);
	std::ofstream outfile("ratios.txt");
	outfile << std::scientific;
	for (auto& R : _R) {
		outfile << R / cgs::GeV << "\t";
		double B = ptr_B10->get_I_R_TOA(R, _phi) + ptr_B11->get_I_R_TOA(R, _phi);
		double C = ptr_C12->get_I_R_TOA(R, _phi);
		double O = ptr_O16->get_I_R_TOA(R, _phi);
		outfile << B / C << "\t";
		outfile << C / O << "\t";
		outfile << "\n";
	}
	outfile.close();
}
