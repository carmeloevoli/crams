#include <iostream>

#include "crams.h"
#include "params.h"

int main() {

	Params par;

	CRAMS crams(par);

	crams.run();

	//crams.test();
	crams.dump();

	return 0;
}



