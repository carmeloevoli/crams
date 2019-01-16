#ifndef INCLUDE_NUCLEILIST_H_
#define INCLUDE_NUCLEILIST_H_

#include <map>

#include "pid.h"

typedef std::map<PID, std::pair<double, double> > NucleiMap;

class NucleiList {
public:
	NucleiList();
	virtual ~NucleiList();
	void add_nucleus(const PID& pid, const double& efficiency, const double& gamma);
	const NucleiMap & get_list() {
		const NucleiMap &ptr = list;
		return ptr;
	}
	size_t get_size() const {
		return list.size();
	}
	NucleiMap list;
};

#endif /* INCLUDE_NUCLEILIST_H_ */
