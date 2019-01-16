#include "output.h"

OutputManager::OutputManager(const LogAxis& T_, const std::vector<Particle>& particles_) :
		T(T_), particles(particles_) {
}

OutputManager::~OutputManager() {
}

void OutputManager::dump_spectra(double R_min, double R_max, size_t R_size) const {
	auto ptr_H1 = find(particles.begin(), particles.end(), Particle(H1, 0));
	auto ptr_C12 = find(particles.begin(), particles.end(), Particle(C12, 0));
	auto ptr_O16 = find(particles.begin(), particles.end(), Particle(O16, 0));
	LogAxis R(R_min, R_max, R_size);
	std::ofstream outfile("primary_spectra.txt");
	outfile << std::scientific;
	double units = 1. / (cgs::GeV * pow2(cgs::meter) * cgs::sec);
	for (size_t i = 0; i < R_size; ++i) {
		double R_ = R.at(i);
		outfile << R_ / cgs::GeV << "\t";
		outfile << ptr_H1->get_I_R_TOA(R_, 0.7 * cgs::GeV) / units << "\t";
		outfile << ptr_C12->get_I_R_TOA(R_, 0.7 * cgs::GeV) / units << "\t";
		outfile << ptr_O16->get_I_R_TOA(R_, 0.7 * cgs::GeV) / units  << "\t";
		outfile << "\n";
	}
	outfile.close();
}

//CRAMS::CRAMS(const Params& par_) :
//		par(par_) {
//	//fill_energy(par.E_min.get(), par.E_max.get(), par.E_size.get());
//	//assert(T.size() == par.E_size.get());
//}
//
//CRAMS::~CRAMS() {
//}
//
//void CRAMS::fill_energy(const double& E_min, const double& E_max, const size_t& E_size) {
//	for (size_t i = 0; i < E_size; ++i) {
//		double ratio = std::pow(E_max / E_min, 1 / (double) (E_size - 1));
//		T.push_back(E_min * std::pow(ratio, (double) i));
//	}
//}
//
///*void CRAMS::fill_rigidity(const double& R_min, const double& R_max, const size_t& R_size) {
// for (size_t i = 0; i < R_size; ++i) {
// double ratio = std::pow(R_max / R_min, 1 / (double) (R_size - 1));
// R.push_back(R_min * std::pow(ratio, (double) i));
// }
// }*/
//
//void CRAMS::fill_particles(const NucleiList& nucleilist) {
////	for (auto nucleus : nucleilist.list) { // TODO pass by reference this https://stackoverflow.com/questions/12851860/accessing-a-map-by-returning-a-pointer-to-it-c
////		Particle particle = Particle(nucleus.first, nucleus.second.first, nucleus.second.second, par);
////		particles.insert(std::make_pair(nucleus.first, particle));
////	}
////	assert(particles.size() == nucleilist.get_size());
//}
//
//void CRAMS::run() {
//	for (auto &p : particles) {
//		std::cout << "# " << p.first << "\n";
//		//p.second.compute_spectrum(T);
//		//p.modulate(T, R);
//	}
//}
//
//void CRAMS::dump() {
//	std::cout << "# R\t";
//	for (auto particle : particles)
//		std::cout << particle.first << "\t";
//	std::cout << "\n";
//	for (size_t i = 0; i < T.size(); ++i) {
//		//std::cout << std::scientific << T.at(i) / GeV << "\t";
//		//for (auto particle : particles)
//		//	std::cout << particle.second.get_I(i) / (1. / GeV / m2 / sec / sr) << "\t";
//		std::cout << "\n";
//	}
//}
