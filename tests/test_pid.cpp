#include <cassert>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "crams/core/pid.h"

using namespace CRAMS;

namespace {

int failures = 0;

void check(bool ok, const char* msg) {
  if (ok) {
    std::cout << "PASS: " << msg << "\n";
  } else {
    std::cerr << "FAIL: " << msg << "\n";
    ++failures;
  }
}

bool approx(double a, double b, double tol = 1e-10) { return std::abs(a - b) < tol; }

// ---------------------------------------------------------------------------
// Construction and getters
// ---------------------------------------------------------------------------

void test_construction() {
  const PID p(6, 12);
  check(p.getZ() == 6, "getZ");
  check(p.getA() == 12, "getA");
  check(p.getId() == 12 * 1000 + 6, "getId = A*1000+Z");
  check(!p.isTertiary(), "isTertiary default false");

  const PID ter(1, 1, true);
  check(ter.isTertiary(), "isTertiary flag");

  // Default-constructed PID is a null state
  const PID null;
  check(null.getZ() == 0, "default Z=0");
  check(null.getA() == 0, "default A=0");
  check(null.getId() == 0, "default id=0");
}

// ---------------------------------------------------------------------------
// getZoverA and getAoverZ
// ---------------------------------------------------------------------------

void test_ZoverA() {
  // Normal nucleus: C12
  check(approx(C12.getZoverA(), 6.0 / 12.0), "C12 Z/A = 0.5");
  check(approx(C12.getAoverZ(), 12.0 / 6.0), "C12 A/Z = 2.0");

  // Proton: Z=1, A=1
  check(approx(H1.getZoverA(), 1.0), "H1 Z/A = 1");
  check(approx(H1.getAoverZ(), 1.0), "H1 A/Z = 1");

  // He4: Z=2, A=4
  check(approx(He4.getZoverA(), 0.5), "He4 Z/A = 0.5");
  check(approx(He4.getAoverZ(), 2.0), "He4 A/Z = 2.0");

  // Default PID (A=0, Z=0): both return 0 without crashing
  const PID null;
  check(approx(null.getZoverA(), 0.0), "null getZoverA = 0");
  check(approx(null.getAoverZ(), 0.0), "null getAoverZ = 0 (no division by zero)");
}

// ---------------------------------------------------------------------------
// Equality and identity
// ---------------------------------------------------------------------------

void test_equality() {
  check(H1 == H1, "H1 == H1");
  check(H1 != H2, "H1 != H2");
  check(H1 != He4, "H1 != He4");
  check(!(H1 == He4), "!(H1 == He4)");

  // Tertiary flag distinguishes otherwise identical PIDs
  check(H1 != H1_ter, "H1 != H1_ter");
  check(H1_ter == H1_ter, "H1_ter == H1_ter");

  // Copies are equal
  const PID copy = C12;
  check(copy == C12, "copy == original");
}

// ---------------------------------------------------------------------------
// Predicate helpers
// ---------------------------------------------------------------------------

void test_predicates() {
  check(H1.isH(), "H1.isH()");
  check(H2.isH(), "H2.isH()");
  check(!He4.isH(), "He4 is not H");
  check(He4.isHe(), "He4.isHe()");
  check(He3.isHe(), "He3.isHe()");
  check(!H1.isHe(), "H1 is not He");
  check(!H1.isTertiary(), "H1 not tertiary");
  check(H1_ter.isTertiary(), "H1_ter is tertiary");
}

// ---------------------------------------------------------------------------
// Ordering: general and special unstable-isotope cases
// ---------------------------------------------------------------------------

void test_ordering() {
  // Heavier nuclei sort after lighter ones (by id = A*1000+Z)
  check(H1 < He4, "H1 < He4");
  check(H1 < C12, "H1 < C12");
  check(C12 < Fe56, "C12 < Fe56");
  check(!(Fe56 < H1), "!(Fe56 < H1)");
  check(!(H1 < H1), "!(H1 < H1) — irreflexive");

  // Same id: tertiary sorts before non-tertiary (so non-tertiary is processed
  // first in the reversed propagation loop, since tertiary sources depend on it)
  check(H1_ter < H1, "H1_ter < H1 (tertiary processed later)");
  check(!(H1 < H1_ter), "!(H1 < H1_ter)");

  // Special cases: unstable isotope sorts after stable counterpart
  // Be10 (unstable) must be processed after Be9 (stable)
  check(!(Be10 < Be9), "Be10 not less than Be9 (processed after)");
  check(Be9 < Be10, "Be9 < Be10");

  // C14 (unstable) after C13 (stable)
  check(!(C14 < C13), "C14 not less than C13");
  check(C13 < C14, "C13 < C14");

  // Cl36 (unstable) after Cl35 (stable)
  check(!(Cl36 < Cl35), "Cl36 not less than Cl35");
  check(Cl35 < Cl36, "Cl35 < Cl36");

  // Mn54 (unstable) after Mn53 (stable)
  check(!(Mn54 < Mn53), "Mn54 not less than Mn53");
  check(Mn53 < Mn54, "Mn53 < Mn54");

  // Verify PID works correctly as a std::map key
  std::map<PID, int> m;
  m[H1] = 1;
  m[He4] = 4;
  m[C12] = 12;
  check(m.size() == 3, "PID usable as map key");
  check(m[C12] == 12, "map lookup by PID");
}

// ---------------------------------------------------------------------------
// toString and operator<<
// ---------------------------------------------------------------------------

void test_string_output() {
  // Both toString() and operator<< should use (A,Z) format
  check(C12.toString() == "(12,6)", "C12.toString() = (12,6)");
  check(H1.toString() == "(1,1)", "H1.toString() = (1,1)");
  check(He4.toString() == "(4,2)", "He4.toString() = (4,2)");

  std::ostringstream oss;
  oss << C12;
  check(oss.str() == "(12,6)", "operator<< C12 = (12,6)");

  std::ostringstream oss2;
  oss2 << H1_ter;
  check(oss2.str() == "(1,1,tertiary)", "operator<< tertiary");

  // toString and operator<< are consistent
  std::ostringstream oss3;
  oss3 << Fe56;
  check(oss3.str() == Fe56.toString(), "toString and operator<< agree");
}

// ---------------------------------------------------------------------------
// Predefined constants spot-check
// ---------------------------------------------------------------------------

void test_constants() {
  check(H1.getZ() == 1 && H1.getA() == 1, "H1: Z=1 A=1");
  check(He4.getZ() == 2 && He4.getA() == 4, "He4: Z=2 A=4");
  check(C12.getZ() == 6 && C12.getA() == 12, "C12: Z=6 A=12");
  check(Fe56.getZ() == 26 && Fe56.getA() == 56, "Fe56: Z=26 A=56");
  check(Be9.getZ() == 4 && Be9.getA() == 9, "Be9: Z=4 A=9");
  check(Be10.getZ() == 4 && Be10.getA() == 10, "Be10: Z=4 A=10");

  // id encoding: A*1000 + Z
  check(Fe56.getId() == 56 * 1000 + 26, "Fe56 id encoding");
}

}  // namespace

int main() {
  test_construction();
  test_ZoverA();
  test_equality();
  test_predicates();
  test_ordering();
  test_string_output();
  test_constants();

  if (failures == 0)
    std::cout << "\nAll " << __FILE__ << " tests passed.\n";
  else
    std::cerr << "\n" << failures << " test(s) failed.\n";

  return failures > 0 ? 1 : 0;
}
