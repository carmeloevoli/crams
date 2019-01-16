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
		return (A > 0) ? fabs((double) Z / (double) A) : 0;
	}
	double get_A_over_Z() const {
		return (A > 0) ? fabs((double) A / (double) Z) : 0;
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
	bool operator>(const PID &other) const {
		return id > other.id;
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

static const PID H1 = PID(1, 1);
static const PID H2 = PID(1, 2);
static const PID Li6 = PID(3, 6);
static const PID Li7 = PID(3, 7);
static const PID Be9 = PID(4, 9);
static const PID Be10 = PID(4, 10);
static const PID B10 = PID(5, 10);
static const PID B11 = PID(5, 11);
static const PID C12 = PID(6, 12);
static const PID C13 = PID(6, 13);
static const PID C14 = PID(6, 14);
static const PID N14 = PID(7, 14);
static const PID N15 = PID(7, 15);
static const PID O16 = PID(8, 16);
static const PID O17 = PID(8, 17);
static const PID O18 = PID(8, 18);
static const PID F19 = PID(9, 19);
static const PID Ne20 = PID(10, 20);
static const PID Ne21 = PID(10, 21);
static const PID Ne22 = PID(10, 22);
static const PID Na22 = PID(11, 22);
static const PID Na23 = PID(11, 23);
static const PID Mg24 = PID(12, 24);
static const PID Mg25 = PID(12, 25);
static const PID Mg26 = PID(12, 26);
static const PID Al26 = PID(13, 26);
static const PID Al27 = PID(13, 27);
static const PID Si28 = PID(14, 28);
static const PID Si29 = PID(14, 29);
static const PID Si30 = PID(14, 30);
static const PID Si32 = PID(14, 32);
static const PID P31 = PID(15, 31);
static const PID P32 = PID(15, 32);
static const PID P33 = PID(15, 33);

//33  	 15 	  2189376 2189376             	  B-      		  P
//32  	 16 	  -1	-1	STABLE                 	          	  S
//33  	 16 	  -1	-1	STABLE                 	          	  S
//34  	 16 	  -1	-1	STABLE                 	          	  S
//35  	 16 	  7560864 7560864              	  B-      		  S
//36  	 16 	  -1	-1	STABLE                 	          	  S
//35  	 17 	  -1	-1	STABLE                 	          	  Cl
//36  	 17 	  9.499E12 9.499E12            	  B-      		  Cl
//37  	 17 	  -1	-1	STABLE                 	          	  Cl
//36  	 18 	  -1	-1	STABLE                 	          	  Ar
//37  	 18 	  3019680.0 3019680.0       	  EC      		  Ar
//38  	 18 	  -1	-1	STABLE                 	          	  Ar
//39  	 18 	  8.489E9 8.489E9             	  B-      		  Ar
//40  	 18 	  -1	-1	STABLE                 	          	  Ar
//42  	 18 	  1.0382E9 1.0382E9            	  B-      		  Ar
//39  	 19 	  -1	-1	STABLE                 	          	  K
//40  	 19 	  3.938E16 3.938E16            	  B-      		  K
//41  	 19 	  -1	-1	STABLE                 	          	  K
//41  	 20 	  3.219E12 3.219E12            	  EC      		  Ca
//42  	 20 	  -1	-1	STABLE                 	          	  Ca
//43  	 20 	  -1	-1	STABLE                 	          	  Ca
//44  	 20 	  -1	-1	STABLE                 	          	  Ca
//45  	 20 	  1.405E7 1.405E7              	  B-      		  Ca
//47  	 20 	  391910.4 391910.4            	  B-      		  Ca
//48  	 20 	  5.996E26 5.996E26            	  BB      		  Ca
//44  	 21 	  210996 210996               	  B-      		  Sc
//45  	 21 	  -1	-1	STABLE                	          	  Sc
//46  	 21 	  7239456.0 7239456.0        	  B-      		  Sc
//47  	 21 	  289370.88 289370.88             B-      		  Sc
//48  	 21 	  157212 157212                   B-      		  Sc
//44  	 22 	  1.893E9 1.893E9                 EC      		  Ti
//46  	 22 	  -1	-1	STABLE                 	          	  Ti
//47  	 22 	  -1	-1	STABLE                 	          	  Ti
//48  	 22 	  -1	-1	STABLE                 	          	  Ti
//49  	 22 	  -1	-1	STABLE                 	          	  Ti
//50  	 22 	  -1	-1	STABLE                 	          	  Ti
//48  	 23 	  1380110.4 1380110.4             EC      		  V
//49  	 23 	  2.851E7 2.851E7                 EC      		  V
//50  	 23 	  4.418E24 4.418E24               EC      		  V
//51  	 23 	  -1	-1	STABLE                 	          	  V
//51  	 24 	  2393496 2393496                 EC      		  Cr
//52  	 24 	  -1	-1	STABLE                 	          	  Cr
//53  	 24 	  -1	-1	STABLE                 	          	  Cr
//54  	 24 	  -1	-1	STABLE                 	          	  Cr
//52  	 25 	  483062.4 483062.4               EC      		  Mn
//53  	 25 	  1.18E14 1.18E14                 EC      		  Mn
//54  	 25 	  2.69817e+07 1.98813e+13         ECB-      		  Mn
//55  	 25 	  -1	-1	STABLE                 	          	  Mn
//54  	 26 	  -1	-1	STABLE                 	          	  Fe
//55  	 26 	  8.61522e+7 8.61522e+7         EC      		  Fe
//56  	 26 	  -1	-1	STABLE                 	          	  Fe
//57  	 26 	  -1	-1	STABLE                 	          	  Fe
//58  	 26 	  -1	-1	STABLE                 	          	  Fe
//59  	 26 	  3844368 3844368                 B-      		  Fe
//60  	 26 	  4.734E13 4.734E13               B-      		  Fe

#endif /* INCLUDE_PID_H_ */
