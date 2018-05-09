#ifndef INCLUDE_PID_H_
#define INCLUDE_PID_H_

#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

class PID {
public:
	PID() {
		set(0, 0);
	}
	PID(const int& Z_, const int& A_) {
		assert(A_ > 0);
		assert(Z_ <= A_);
		set(Z_, A_);
	}
	virtual ~PID() {
	}
	void set(const int& Z_, const int& A_) {
		Z = Z_;
		A = A_;
		id = A * 1000 + Z;
	}
	int get_Z() const {
		return Z;
	}
	int get_A() const {
		return A;
	}
	double get_Z_over_A() const {
		return (A > 0) ? fabs((double)Z / (double)A) : 0;
	}
	double get_A_over_Z() const {
		return (A > 0) ? fabs((double)A / (double)Z) : 0;
	}
	int get_id() const {
		return id;
	}
	bool operator==(const PID &other) const {
		return id == other.id;
	}
	bool operator!=(const PID &other) const {
		return id != other.id;
	}
	bool operator<(const PID &other) const {
		return id < other.id;
	}
	bool with_Z(const int& Z_) const {
		return Z == Z_;
	}
	bool with_A(const int& A_) const {
		return A == A_;
	}
	bool is_lepton() const {
		return A == 0;
	}
	bool is_electron() const {
		return A == 0 && Z == -1;
	}
	bool is_positron() const {
		return A == 0 && Z == 1;
	}
	bool is_H() const {
		return (Z == 1);
	}
	bool is_He() const {
		return (Z == 2);
	}
	friend std::ostream& operator<<(std::ostream& stream, const PID& pid) {
		stream << "(" << pid.get_A() << "," << pid.get_Z() << ")";
		return stream;
	}
	std::string to_string() const {
		std::string ss;
		ss = "(" + std::to_string(Z) + "," + std::to_string(A) + ")";
		return ss;
	}
protected:
	int Z;
	int A;
	int id;
};

#endif /* INCLUDE_PID_H_ */
