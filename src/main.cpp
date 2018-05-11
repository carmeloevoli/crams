#include <iostream>

#include "chi2.h"
#include "crams.h"
#include "nucleilist.h"
#include "params.h"
#include "pid.h"

int main() {

	Params par;
	par.potential.set(0.1 * GeV);
	par.v_A.set(1 * km / sec);
	par.D_0.set(3e29 * cm2 / sec);
	par.R_b.set(100 * GeV);
	par.delta_s.set(0.1);
	par.delta_low.set(0.5);
	par.print();

	NucleiList nucleilist;
	//nucleilist.add_nucleus(PID(1, 1), 0.05, 4.2);
	//nucleilist.add_nucleus(PID(2, 4), 0.1, 4.2);
	nucleilist.add_nucleus(PID(6, 12), 0.01, 4.3);
	nucleilist.add_nucleus(PID(8, 16), 0.01, 4.3);

	Chi2 Carbon("C_AMS02_rigidity.txt");
	Chi2 Oxygen("O_AMS02_rigidity.txt");

	CRAMS crams(par);
	crams.fill_particles(nucleilist);
	crams.run();
	crams.dump();

	PID C12 = PID(6, 12);
	auto T = crams.get_T();

	Carbon.calculate_chi2(C12, T, crams.get_spectrum(C12), par.potential.get());

	std::cout << "\n" << Carbon.getChi2() << "\n";

	return 0;
}



