#include <iostream>
#include <string>

#include "params.h"

ParticleList::~ParticleList() {
	list.clear();
	std::cout << "released memory from ParticleList\n";
}

bool ParticleList::insert(const PID& key, const double& value) {
	auto res = list.insert(std::make_pair(key, value));
	if (!res.second) {
		std::cout << "particle " << key << " already exists " << " with abundance "
				<< (res.first)->second << std::endl;
	} else {
		std::cout << "created particle " << key << " with abundance " << value << std::endl;
	}
	return res.second;
}

void ParticleList::set_abundance(const PID& key, const double& value) {
	auto it = list.find(key);
	if (it != list.end()) {
		it->second = value;
		std::cout << "PID : " << key << " abundance modified to " << value << ".\n";
	} else
		std::cout << "PID : " << key << " not found in particle list.\n";
}

void ParticleList::print() {
	std::cout << "Particle list contains " << list.size() << " nuclei.\n";
	for (auto& particle : list)
		std::cout << particle.first << "\n";
}

Params::~Params() {
	std::cout << "released memory from Params\n";
}

void Params::print() {
	std::cout << "H     : " << _H / cgs::kpc << " kpc\n";
	std::cout << "E_min : " << _T_min / cgs::GeV << " GeV\n";
	std::cout << "E_max : " << _T_max / cgs::GeV << " GeV\n";
}
