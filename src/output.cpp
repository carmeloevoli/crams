#include "output.h"
#include "utilities.h"

OutputManager::OutputManager(const Particles& particles, const double phi) :
		_particles(particles), _phi(phi) {
}

OutputManager::~OutputManager() {
}

ptr_Particle OutputManager::find_ptr(const PID& pid) {
	auto it = find(_particles.begin(), _particles.end(), Particle(pid, 0));
	bool isPresent = !(it == _particles.end());
	return {isPresent, it};
}

double OutputManager::B(const double& R) const {
	double value = (ptr_B10.isPresent) ? ptr_B10.it->I_R_TOA(R, _phi) : 0;
	value += (ptr_B11.isPresent) ? ptr_B11.it->I_R_TOA(R, _phi) : 0;
	return value;
}

double OutputManager::C(const double& R) const {
	double value = (ptr_C12.isPresent) ? ptr_C12.it->I_R_TOA(R, _phi) : 0.;
	value += (ptr_C13.isPresent) ? ptr_C13.it->I_R_TOA(R, _phi) : 0.;
	value += (ptr_C14.isPresent) ? ptr_C14.it->I_R_TOA(R, _phi) : 0.;
	return value;
}

double OutputManager::N(const double& R) const {
	double value = (ptr_N14.isPresent) ? ptr_N14.it->I_R_TOA(R, _phi) : 0.;
	value += (ptr_N15.isPresent) ? ptr_N15.it->I_R_TOA(R, _phi) : 0.;
	return value;
}

double OutputManager::O(const double& R) const {
	double value = (ptr_O16.isPresent) ? ptr_O16.it->I_R_TOA(R, _phi) : 0;
	value += (ptr_O17.isPresent) ? ptr_O17.it->I_R_TOA(R, _phi) : 0.;
	value += (ptr_O18.isPresent) ? ptr_O18.it->I_R_TOA(R, _phi) : 0.;
	return value;
}

double OutputManager::C12_C13(const double& R) const {
	double C12 = (ptr_C12.isPresent) ? ptr_C12.it->I_R_TOA(R, _phi) : 0.;
	double C13 = (ptr_C13.isPresent) ? ptr_C13.it->I_R_TOA(R, _phi) : 0.;
	return C12 / C13;
}

void OutputManager::dump_spectra(double R_min, double R_max, size_t R_size) const {
	auto _R = LogAxis(R_min, R_max, R_size);
	std::ofstream outfile("spectra.txt");
	outfile << std::scientific;
	double units = 1. / (cgs::GeV * pow2(cgs::meter) * cgs::sec);
	for (auto& R : _R) {
		outfile << R / cgs::GeV << "\t";
		outfile << 0 / units << "\t";
		outfile << B(R) / units << "\t";
		outfile << C(R) / units << "\t";
		outfile << N(R) / units << "\t";
		outfile << O(R) / units << "\t";
		outfile << "\n";
	}
	outfile.close();
}

void OutputManager::dump_ratio(double R_min, double R_max, size_t R_size) const {
	auto _R = LogAxis(R_min, R_max, R_size);
	std::ofstream outfile("ratios.txt");
	outfile << std::scientific;
	for (auto& R : _R) {
		outfile << R / cgs::GeV << "\t";
		outfile << B(R) / C(R) << "\t";
		outfile << C(R) / O(R) << "\t";
		outfile << B(R) / O(R) << "\t";
		outfile << "\n";
	}
	outfile.close();
}

void OutputManager::dump_heavy_spectra(double R_min, double R_max, size_t R_size) const {
//	auto ptr_Ne20 = find(_particles.begin(), _particles.end(), Particle(Ne20, 0));
//	auto ptr_Mg24 = find(_particles.begin(), _particles.end(), Particle(Mg24, 0));
//	auto ptr_Si28 = find(_particles.begin(), _particles.end(), Particle(Si28, 0));
//	auto ptr_Fe56 = find(_particles.begin(), _particles.end(), Particle(Fe56, 0));
//	auto _R = LogAxis(R_min, R_max, R_size);
//	std::ofstream outfile("heavy_spectra.txt");
//	outfile << std::scientific;
//	double units = 1. / (cgs::GeV * pow2(cgs::meter) * cgs::sec);
//	for (auto& R : _R) {
//		outfile << R / cgs::GeV << "\t";
//		outfile << ptr_Ne20->I_R_TOA(R, _phi) / units << "\t";
//		outfile << ptr_Mg24->I_R_TOA(R, _phi) / units << "\t";
//		outfile << ptr_Si28->I_R_TOA(R, _phi) / units << "\t";
//		outfile << ptr_Fe56->I_R_TOA(R, _phi) / units << "\t";
//		outfile << "\n";
//	}
//	outfile.close();
}
