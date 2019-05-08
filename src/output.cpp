#include "output.h"
#include "utilities.h"

OutputManager::OutputManager(const Particles& particles, double phi, size_t id) :
		_particles(particles), _phi(phi) {
	spectra_filename = "spectra_" + std::to_string(id) + ".txt";
	ratios_filename = "ratios_" + std::to_string(id) + ".txt";
}

OutputManager::~OutputManager() {
}

ptr_Particle OutputManager::find_ptr(const PID& pid) {
	auto it = find(_particles.begin(), _particles.end(), Particle(pid, 0));
	bool isPresent = !(it == _particles.end());
	return {isPresent, it};
}

double OutputManager::H(const double& R) const {
	double value = (ptr_H1.isPresent) ? ptr_H1.it->I_R_TOA(R, _phi) : 0;
	value += (ptr_H1_ter.isPresent) ? ptr_H1_ter.it->I_R_TOA(R, _phi) : 0;
	value += (ptr_H2.isPresent) ? ptr_H2.it->I_R_TOA(R, _phi) : 0;
	return value;
}

double OutputManager::He(const double& R) const {
	double value = (ptr_He4.isPresent) ? ptr_He4.it->I_R_TOA(R, _phi) : 0;
	value += (ptr_He3.isPresent) ? ptr_He3.it->I_R_TOA(R, _phi) : 0;
	return value;
}

double OutputManager::Li(const double& R) const {
	double value = (ptr_Li6.isPresent) ? ptr_Li6.it->I_R_TOA(R, _phi) : 0;
	value += (ptr_Li7.isPresent) ? ptr_Li7.it->I_R_TOA(R, _phi) : 0;
	return value;
}

double OutputManager::Be(const double& R) const {
	double value = (ptr_Be9.isPresent) ? ptr_Be9.it->I_R_TOA(R, _phi) : 0;
	value += (ptr_Be10.isPresent) ? ptr_Be10.it->I_R_TOA(R, _phi) : 0;
	value += (ptr_Be7.isPresent) ? ptr_Be7.it->I_R_TOA(R, _phi) : 0;
	return value;
}

double OutputManager::B(const double& R) const {
	double value = (ptr_B10.isPresent) ? ptr_B10.it->I_R_TOA(R, _phi) : 0;
	value += (ptr_B11.isPresent) ? ptr_B11.it->I_R_TOA(R, _phi) : 0;
	return value;
}

double OutputManager::C(const double& R) const { // TODO TOA -> LIS
	double value = (ptr_C12.isPresent) ? ptr_C12.it->I_R_LIS(R) : 0.;
	value += (ptr_C13.isPresent) ? ptr_C13.it->I_R_LIS(R) : 0.;
	value += (ptr_C14.isPresent) ? ptr_C14.it->I_R_LIS(R) : 0.;
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

void OutputManager::dump_spectra(double R_min, double R_max, size_t R_size) const {
	auto _R = LogAxis(R_min, R_max, R_size);
	std::ofstream outfile(spectra_filename.c_str());
	outfile << std::scientific;
	double units = 1. / (cgs::GeV * pow2(cgs::meter) * cgs::sec);
	for (auto& R : _R) {
		outfile << R / cgs::GeV << "\t";
		outfile << H(R) / units << "\t";
		outfile << He(R) / units << "\t";
		outfile << Li(R) / units << "\t";
		outfile << Be(R) / units << "\t";
		outfile << B(R) / units << "\t";
		outfile << C(R) / units << "\t";
		outfile << N(R) / units << "\t";
		outfile << O(R) / units << "\t";
		outfile << "\n";
	}
	outfile.close();
}

void OutputManager::dump_ratios(double R_min, double R_max, size_t R_size) const {
	auto _R = LogAxis(R_min, R_max, R_size);
	std::ofstream outfile(ratios_filename.c_str());
	outfile << std::scientific;
	for (auto& R : _R) {
		outfile << R / cgs::GeV << "\t";
		outfile << Li(R) / B(R) << "\t";
		outfile << Be(R) / B(R) << "\t";
		outfile << Li(R) / C(R) << "\t";
		outfile << Be(R) / C(R) << "\t";
		outfile << B(R) / C(R) << "\t";
		outfile << C(R) / O(R) << "\t";
		outfile << B(R) / O(R) << "\t";
		outfile << He(R) / O(R) << "\t";
		outfile << "\n";
	}
	outfile.close();
}

void OutputManager::dump_heavy_spectra(double R_min, double R_max, size_t R_size) const {
}
