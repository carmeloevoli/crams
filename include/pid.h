#ifndef INCLUDE_PID_H_
#define INCLUDE_PID_H_

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

//#include "spdlog/fmt/ostr.h"  // must be included
//#include "spdlog/spdlog.h"    // must be included

namespace CRAMS {

class PID {
 public:
  PID() { set(0, 0, false); }

  PID(const int& Z, const int& A, const bool& isTertiary = false) {
    assert(A > 0);
    assert(Z <= A);
    set(Z, A, isTertiary);
  }

  virtual ~PID() {}

  void set(const int& Z, const int& A, const bool& isTertiary) {
    m_Z = Z;
    m_A = A;
    m_id = A * 1000 + Z;
    m_isTertiary = isTertiary;
  }

  int getZ() const { return m_Z; }
  int getA() const { return m_A; }
  double getZoverA() const { return (m_A > 0) ? fabs((double)m_Z / (double)m_A) : 0; }
  double getAoverZ() const { return (m_A > 0) ? fabs((double)m_A / (double)m_Z) : 0; }
  int getId() const { return m_id; }

  bool operator==(const PID& other) const { return m_id == other.m_id && m_isTertiary == other.m_isTertiary; }
  bool operator!=(const PID& other) const { return m_id != other.m_id || m_isTertiary != other.m_isTertiary; }

  bool operator<(const PID& other) const {
    // Be10
    if (m_id == 10005 && other.m_id == 10004) return true;
    if (m_id == 10004 && other.m_id == 10005) return false;
    // C14
    if (m_id == 14007 && other.m_id == 14006) return true;
    if (m_id == 14006 && other.m_id == 14007) return false;
    // Cl36
    if (m_id == 54026 && other.m_id == 54025) return true;
    if (m_id == 54025 && other.m_id == 54026) return false;
    // Mn54
    if (m_id == 36018 && other.m_id == 36017) return true;
    if (m_id == 36017 && other.m_id == 36018) return false;
    // others
    if (m_id != other.m_id)
      return m_id < other.m_id;
    else
      return m_isTertiary;
  }

  bool isH() const { return (m_Z == 1); }
  bool isHe() const { return (m_Z == 2); }
  bool isTertiary() const { return m_isTertiary; }

  friend std::ostream& operator<<(std::ostream& stream, const PID& pid) {
    stream << "(" << pid.getA() << "," << pid.getZ() << ")";
    return stream;
  }

  // template <typename OStream>
  // friend OStream& operator<<(OStream& os, const PID& c) {
  //   return os << "pid";
  // }

  std::string toString() const {
    std::string ss;
    ss = "(" + std::to_string(m_Z) + "," + std::to_string(m_A) + ")";
    return ss;
  }

 protected:
  int m_Z;
  int m_A;
  int m_id;
  bool m_isTertiary;
};  // namespace CRAMS

typedef std::pair<PID, PID> Channel;

static const PID H1_ter = PID(1, 1, true);
static const PID H1 = PID(1, 1);
static const PID H2 = PID(1, 2);
static const PID He3 = PID(2, 3);
static const PID He4 = PID(2, 4);
static const PID Li6 = PID(3, 6);
static const PID Li7 = PID(3, 7);
static const PID Be7 = PID(4, 7);
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
static const PID S32 = PID(16, 32);
static const PID S33 = PID(16, 33);
static const PID S34 = PID(16, 34);
static const PID S36 = PID(16, 36);
static const PID Cl35 = PID(17, 35);
static const PID Cl36 = PID(17, 36);
static const PID Cl37 = PID(17, 37);
static const PID Ar36 = PID(18, 36);
static const PID Ar37 = PID(18, 37);
static const PID Ar38 = PID(18, 38);
static const PID Ar40 = PID(18, 40);
static const PID K39 = PID(19, 39);
static const PID K40 = PID(19, 40);
static const PID K41 = PID(19, 41);
static const PID Ca40 = PID(20, 40);
static const PID Ca41 = PID(20, 41);
static const PID Ca42 = PID(20, 42);
static const PID Ca43 = PID(20, 43);
static const PID Ca44 = PID(20, 44);
static const PID Ca46 = PID(20, 46);
static const PID Ca48 = PID(20, 48);
static const PID Sc45 = PID(21, 45);
static const PID Ti44 = PID(22, 44);
static const PID Ti46 = PID(22, 46);
static const PID Ti47 = PID(22, 47);
static const PID Ti48 = PID(22, 48);
static const PID Ti49 = PID(22, 49);
static const PID Ti50 = PID(22, 50);
static const PID V49 = PID(23, 49);
static const PID V50 = PID(23, 50);
static const PID V51 = PID(23, 51);
static const PID Cr48 = PID(24, 48);
static const PID Cr50 = PID(24, 50);
static const PID Cr51 = PID(24, 51);
static const PID Cr52 = PID(24, 52);
static const PID Cr53 = PID(24, 53);
static const PID Cr54 = PID(24, 54);
static const PID Mn53 = PID(25, 53);
static const PID Mn54 = PID(25, 54);
static const PID Mn55 = PID(25, 55);
static const PID Fe54 = PID(26, 54);
static const PID Fe55 = PID(26, 55);
static const PID Fe56 = PID(26, 56);
static const PID Fe57 = PID(26, 57);
static const PID Fe58 = PID(26, 58);
static const PID Fe60 = PID(26, 60);
static const PID Co57 = PID(27, 57);
static const PID Co59 = PID(27, 59);
static const PID Ni56 = PID(28, 56);
static const PID Ni58 = PID(28, 58);
static const PID Ni59 = PID(28, 59);
static const PID Ni60 = PID(28, 60);
static const PID Ni61 = PID(28, 61);
static const PID Ni62 = PID(28, 62);
static const PID Ni64 = PID(28, 64);

}  // namespace CRAMS

#endif  // INCLUDE_PID_H_
