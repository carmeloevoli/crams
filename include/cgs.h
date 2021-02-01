#ifndef INCLUDE_CGS_H_
#define INCLUDE_CGS_H_

#include <cmath>

namespace CRAMS {
namespace CGS {

// CGS UNITS
static constexpr double second = 1;
static constexpr double centimeter = 1.;
static constexpr double gram = 1;
static constexpr double kelvin = 1;
static constexpr double sr = 1;

// TIME UNITS
static constexpr double year = 3.154e+7 * second;
static constexpr double kiloyear = 1e3 * year;
static constexpr double Megayear = 1e6 * year;
static constexpr double Gigayear = 1e9 * year;

// LENGTH UNITS
static constexpr double meter = 1e2 * centimeter;
static constexpr double kilometer = 1e3 * meter;
static constexpr double parsec = 3.086e16 * meter;
static constexpr double kiloparsec = 1e3 * parsec;
static constexpr double fm = 1e-13 * centimeter;

// MASS UNITS
static constexpr double mgram = 1e-3 * gram;
static constexpr double kilogram = 1e3 * gram;

// ENERGY UNITS
static constexpr double erg = gram * centimeter * centimeter / second;
static constexpr double joule = 1e7 * erg;
static constexpr double electronvolt = 1.60217657e-19 * joule;
static constexpr double kiloelectronvolt = 1e3 * electronvolt;
static constexpr double megaelectronvolt = 1e6 * electronvolt;
static constexpr double gigaelectronvolt = 1e9 * electronvolt;
static constexpr double teraelectronvolt = 1e12 * electronvolt;
static constexpr double petaelectronvolt = 1e15 * electronvolt;

// ABBREVIATION
static constexpr double sec = second;
static constexpr double km = kilometer;
static constexpr double kyr = kiloyear;
static constexpr double Myr = Megayear;
static constexpr double kpc = kiloparsec;
static constexpr double eV = electronvolt;
static constexpr double keV = kiloelectronvolt;
static constexpr double MeV = megaelectronvolt;
static constexpr double GeV = gigaelectronvolt;
static constexpr double TeV = teraelectronvolt;
static constexpr double PeV = petaelectronvolt;
static constexpr double cm = centimeter;
static constexpr double cm2 = cm * cm;
static constexpr double cm3 = cm * cm * cm;
static constexpr double m2 = meter * meter;

// PHYSICAL CONSTANTS
static constexpr double cLight = 2.99792458e10 * centimeter / second;
static constexpr double cSquared = cLight * cLight;
static constexpr double protonMass = 1.67262158e-24 * gram;
static constexpr double protonMassC = protonMass * cLight;
static constexpr double protonMassC2 = protonMass * cSquared;
static constexpr double neutronMass = 1.67492735e-24 * gram;
static constexpr double neutronMassC2 = neutronMass * cSquared;
static constexpr double electronMass = 9.10938291e-28 * gram;
static constexpr double electronMassC2 = electronMass * cSquared;
static constexpr double sunMass = 1.989e33 * gram;
static constexpr double hPlanck = 6.62607015e-34 * joule * second;
static constexpr double kBoltzmann = 1.3806488e-23 * joule / kelvin;
static constexpr double electronRadius = 2.8179403227e-15 * meter;
static constexpr double IsH = 19 * eV;   // H  eff. ioniz. potential
static constexpr double IsHe = 44 * eV;  // He eff. ioniz. potential
static constexpr double barn = 1e-24 * cm2;
static constexpr double mbarn = 1e-3 * barn;

// MOMENTUM UNITS
static constexpr double eV_c = electronvolt / cLight;
static constexpr double keV_c = 1e3 * eV_c;
static constexpr double MeV_c = 1e6 * eV_c;
static constexpr double GeV_c = 1e9 * eV_c;
static constexpr double TeV_c = 1e12 * eV_c;
static constexpr double PeV_c = 1e15 * eV_c;

// SPECIFIC CODE COSTANTS
static constexpr double E_SN = 1e51 * erg;
static constexpr double snRate = 1. / (30. * year);
static constexpr double galaxySize = 10. * kpc;
static constexpr double f_He = 0.1;
static constexpr double K_He = 2.51984209979;  // 4**2/3
static constexpr double meanISMmass = protonMass * (1. + 4 * f_He) / (1. + f_He);
static constexpr double inelasticity = 0.5;

}  // namespace CGS
}  // namespace CRAMS

#endif  // INCLUDE_CGS_H_
