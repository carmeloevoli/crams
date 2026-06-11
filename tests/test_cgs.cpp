#include <cmath>
#include <iostream>

#include "crams/core/cgs.h"

using namespace CRAMS::CGS;

// ---------------------------------------------------------------------------
// Compile-time checks (exact definitions and ordering invariants)
// ---------------------------------------------------------------------------

// Derived quantities match their definitions
static_assert(cSquared == cLight * cLight, "cSquared");
static_assert(protonMassC2 == protonMass * cSquared, "protonMassC2");
static_assert(neutronMassC2 == neutronMass * cSquared, "neutronMassC2");
static_assert(electronMassC2 == electronMass * cSquared, "electronMassC2");
static_assert(m2 == meter * meter, "m2");
static_assert(cm2 == cm * cm, "cm2");
static_assert(cm3 == cm * cm * cm, "cm3");
static_assert(barn == 1e-24 * cm2, "barn");
static_assert(mbarn == 1e-3 * barn, "mbarn");

// Abbreviations match their full names
static_assert(sec == second, "sec alias");
static_assert(km == kilometer, "km alias");
static_assert(kpc == kiloparsec, "kpc alias");
static_assert(eV == electronvolt, "eV alias");
static_assert(keV == kiloelectronvolt, "keV alias");
static_assert(MeV == megaelectronvolt, "MeV alias");
static_assert(GeV == gigaelectronvolt, "GeV alias");
static_assert(TeV == teraelectronvolt, "TeV alias");
static_assert(PeV == petaelectronvolt, "PeV alias");

// Ordering invariants
static_assert(neutronMassC2 > protonMassC2, "neutron heavier than proton");
static_assert(protonMass > electronMass, "proton heavier than electron");
static_assert(protonMassC2 > electronMassC2, "proton rest energy > electron rest energy");
static_assert(GeV > MeV, "GeV > MeV");
static_assert(MeV > keV, "MeV > keV");
static_assert(TeV > GeV, "TeV > GeV");
static_assert(kpc > parsec, "kpc > parsec");
static_assert(Gigayear > Megayear, "Gyr > Myr");
static_assert(Megayear > kiloyear, "Myr > kyr");

// ---------------------------------------------------------------------------
// Runtime checks (approximate physical values against known literature)
// ---------------------------------------------------------------------------

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

bool approx(double a, double b, double tol = 1e-4) { return std::abs(a / b - 1.0) < tol; }

void test_particle_masses() {
  check(approx(protonMassC2 / MeV, 938.272), "proton rest mass = 938.272 MeV");
  check(approx(neutronMassC2 / MeV, 939.565), "neutron rest mass = 939.565 MeV");
  check(approx(electronMassC2 / MeV, 0.510999), "electron rest mass = 0.511 MeV");
  // neutron-proton mass difference ~ 1.293 MeV
  check(approx((neutronMassC2 - protonMassC2) / MeV, 1.293, 1e-3), "n-p mass difference = 1.293 MeV");
}

void test_speed_of_light() { check(approx(cLight, 2.99792458e10 * cm / sec), "c = 2.998e10 cm/s"); }

void test_unit_chains() {
  check(approx(GeV / MeV, 1e3), "GeV / MeV = 1000");
  check(approx(MeV / keV, 1e3), "MeV / keV = 1000");
  check(approx(TeV / GeV, 1e3), "TeV / GeV = 1000");
  check(approx(PeV / TeV, 1e3), "PeV / TeV = 1000");
  check(approx(kilometer / meter, 1e3), "km / m = 1000");
  check(approx(meter / centimeter, 1e2), "m / cm = 100");
  check(approx(kpc / parsec, 1e3), "kpc / parsec = 1000");
  check(approx(kiloyear / year, 1e3), "kyr / yr = 1000");
  check(approx(Megayear / kiloyear, 1e3), "Myr / kyr = 1000");
}

void test_length_units() {
  // 1 parsec ~ 3.086e18 cm
  check(approx(parsec / (3.086e18 * cm), 1.0, 1e-3), "parsec ~ 3.086e18 cm");
  // 1 fm = 1e-13 cm
  check(approx(fm, 1e-13 * cm), "fm = 1e-13 cm");
}

void test_energy_units() {
  // 1 eV in erg
  check(approx(eV / erg, 1.602176634e-12), "eV = 1.602e-12 erg");
  // 1 joule = 1e7 erg
  check(approx(joule / erg, 1e7), "joule = 1e7 erg");
}

void test_cross_sections() {
  check(approx(barn, 1e-24 * cm2), "barn = 1e-24 cm2");
  check(approx(mbarn / barn, 1e-3), "mbarn = 1e-3 barn");
}

}  // namespace

int main() {
  test_particle_masses();
  test_speed_of_light();
  test_unit_chains();
  test_length_units();
  test_energy_units();
  test_cross_sections();

  if (failures == 0)
    std::cout << "\nAll " << __FILE__ << " tests passed.\n";
  else
    std::cerr << "\n" << failures << " test(s) failed.\n";

  return failures > 0 ? 1 : 0;
}
