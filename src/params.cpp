#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "params.h"

ParticleList::~ParticleList() {
	_list.clear();
	std::cout << "released memory from ParticleList\n";
}

bool ParticleList::insert(const PID& key, const double& value) {
	auto res = _list.insert(std::make_pair(key, value));
	if (!res.second) {
		std::cout << "particle " << key << " already exists " << " with abundance "
				<< (res.first)->second << std::endl;
	} else {
		std::cout << "created particle " << key << " with abundance " << value << std::endl;
	}
	return res.second;
}

void ParticleList::set_abundance(const PID& key, const double& value) {
	auto it = _list.find(key);
	if (it != _list.end()) {
		it->second = value;
		std::cout << "PID : " << key << " abundance modified to " << value << ".\n";
	} else
		std::cout << "PID : " << key << " not found in particle list.\n";
}

void ParticleList::set_from_file(const std::string& filename) {
	std::ifstream infile(filename.c_str());
	std::string line;
	while (std::getline(infile, line)) {
		std::istringstream iss(line);
		std::string key;
		double value;
		if (!(iss >> key >> value)) {
			break;
		} // error
		if (key == "q_H")
			set_abundance(H1, value);
		else if (key == "q_C")
			set_abundance(C12, value);
		else if (key == "q_N")
			set_abundance(N14, value);
		else if (key == "q_O")
			set_abundance(O16, value);
	}
}

void ParticleList::print() {
	std::cout << "Particle list contains " << _list.size() << " nuclei.\n";
	for (auto& particle : _list)
		std::cout << particle.first << "\n";
}

Params::~Params() {
	std::cout << "released memory from Params\n";
}

void Params::set_params(const std::string& key, const double& value) {
	if (key == "D_0")
		_D_0 = value * 1e28 * cgs::cm2 / cgs::sec;
	else if (key == "H")
		_H = value * cgs::kpc;
	else if (key == "delta")
		_delta = value;
	else if (key == "ddelta")
		_ddelta = value;
	else if (key == "R_b")
		_R_b = value * cgs::GeV;
	else if (key == "v_A")
		_v_A = value * cgs::km / cgs::sec;
	else if (key == "H_slope")
		_H_slope = value;
	else if (key == "slope")
		_nuclei_slope = value;
	else if (key == "phi")
		_modulation_potential = value * cgs::GeV;
	else if (key == "id")
		_id = (int) value;
}

void Params::set_from_file(const std::string& filename) {
	std::ifstream infile(filename.c_str());
	std::string line;
	while (std::getline(infile, line)) {
		std::istringstream iss(line);
		std::string key;
		double value;
		if (!(iss >> key >> value)) {
			break;
		} // error
		set_params(key, value);
	}
}

void Params::print() {
	std::cout << "H      : " << _H / cgs::kpc << " kpc\n";
	std::cout << "mu     : " << _mu / (cgs::mgram / cgs::cm2) << " mg/cm2\n";
	std::cout << "vA     : " << _v_A / (cgs::km / cgs::sec) << " km/s\n";
	std::cout << "D0     : " << _D_0 / (1e28 * cgs::cm2 / cgs::sec) << " x 1e28 cm2/s\n";
	std::cout << "delta  : " << _delta << "\n";
	std::cout << "ddelta : " << _ddelta << "\n";
	std::cout << "R_b    : " << _R_b / cgs::GeV << " GeV\n";
	std::cout << "s      : " << _smoothness << "\n";
	std::cout << "H_slope: " << _H_slope << "\n";
	std::cout << "slope  : " << _nuclei_slope << "\n";
	std::cout << "phi    : " << _modulation_potential / cgs::GeV << " GeV\n";
	std::cout << "E_min  : " << _T_min / cgs::GeV << " GeV\n";
	std::cout << "E_max  : " << _T_max / cgs::GeV << " GeV\n";
	std::cout << "E_size : " << _T_size << "\n";
}
