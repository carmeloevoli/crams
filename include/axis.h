#ifndef INCLUDE_AXIS_H_
#define INCLUDE_AXIS_H_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

template<class T>
class TAxis {

protected:
	T min;
	T max;
	size_t size;
	std::vector<T> axis;

public:
	TAxis() :
			min(0), max(0), size(0) {
	}

	virtual ~TAxis() {
		axis.clear();
	}

	/* Accessor / Mutator */
	T &get(size_t i) {
		return axis[i];
	}

	T &get(const size_t& i) {
		return axis[i];
	}

	/* Accessor */
	const T &get(size_t i) const {
		return axis[i];
	}

	const T &get(const size_t& i) const {
		return axis[i];
	}

	T at(size_t i) const {
		return axis[i];
	}

	T get_max() const {
		return max;
	}

	T get_min() const {
		return min;
	}

	size_t get_size() const {
		return size;
	}

	const std::vector<T>& get_axis() const {
		return axis;
	}

	void show_axis(const std::string& name, const T& units) const {
		auto min = *min_element(axis.begin(), axis.end());
		auto max = *max_element(axis.begin(), axis.end());
		std::cout << name << " axis in : ";
		std::cout << min / units << " ... " << max / units;
		std::cout << " with " << axis.size() << " elements.\n";
	}

	virtual void set(const T& min, const T& max, const size_t& size) {
		this->min = min;
		this->max = max;
		this->size = size;
	}
};

class LinAxis: public TAxis<double> {
public:
	LinAxis(const double& min, const double& max, const size_t& size) {
		set(min, max, size);
		const double dx = (max - min) / (double) (size - 1);
		for (size_t i = 0; i < size; ++i) {
			double value = min + dx * i;
			axis.push_back(value);
		}
	}
};

class LogAxis: public TAxis<double> {
public:
	LogAxis(const double& min, const double& max, const size_t& size) {
		set(min, max, size);
		const double delta_log = std::exp(std::log(max / min) / (size - 1));
		for (size_t i = 0; i < size; ++i) {
			double value = std::exp(std::log(min) + (double) i * std::log(delta_log));
			axis.push_back(value);
		}
	}

	inline double get_ratio() const {
		return axis[1] / axis[0];
	}
};

#endif /* INCLUDE_AXIS_H_ */
