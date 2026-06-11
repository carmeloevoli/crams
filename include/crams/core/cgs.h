#ifndef CRAMS_CORE_CGS_H_
#define CRAMS_CORE_CGS_H_

#include <cmath>

namespace CRAMS {
namespace CGS {

// CGS UNITS
constexpr double second = 1;
constexpr double centimeter = 1.;
constexpr double gram = 1;
constexpr double kelvin = 1;
constexpr double sr = 1;

// TIME UNITS
constexpr double year = 3.15576e7 * second;  // Julian year: 365.25 * 86400 s
constexpr double kiloyear = 1e3 * year;
constexpr double Megayear = 1e6 * year;
constexpr double Gigayear = 1e9 * year;

// LENGTH UNITS
constexpr double meter = 1e2 * centimeter;
constexpr double kilometer = 1e3 * meter;
constexpr double parsec = 3.085677581e16 * meter;  // IAU 2012
constexpr double kiloparsec = 1e3 * parsec;
constexpr double fm = 1e-13 * centimeter;

// MASS UNITS
constexpr double mgram = 1e-3 * gram;
constexpr double kilogram = 1e3 * gram;

// ENERGY UNITS
constexpr double erg = gram * centimeter * centimeter / (second * second);
constexpr double joule = 1e7 * erg;
constexpr double electronvolt = 1.602176634e-19 * joule;  // CODATA 2018 exact
constexpr double kiloelectronvolt = 1e3 * electronvolt;
constexpr double megaelectronvolt = 1e6 * electronvolt;
constexpr double gigaelectronvolt = 1e9 * electronvolt;
constexpr double teraelectronvolt = 1e12 * electronvolt;
constexpr double petaelectronvolt = 1e15 * electronvolt;

// ABBREVIATIONS
constexpr double sec = second;
constexpr double km = kilometer;
constexpr double kyr = kiloyear;
constexpr double Myr = Megayear;
constexpr double kpc = kiloparsec;
constexpr double eV = electronvolt;
constexpr double keV = kiloelectronvolt;
constexpr double MeV = megaelectronvolt;
constexpr double GeV = gigaelectronvolt;
constexpr double TeV = teraelectronvolt;
constexpr double PeV = petaelectronvolt;
constexpr double cm = centimeter;
constexpr double cm2 = cm * cm;
constexpr double cm3 = cm * cm * cm;
constexpr double m2 = meter * meter;

// PHYSICAL CONSTANTS (CODATA 2018)
constexpr double cLight = 2.99792458e10 * centimeter / second;  // exact
constexpr double cSquared = cLight * cLight;
constexpr double protonMass = 1.67262192369e-24 * gram;
constexpr double protonMassC = protonMass * cLight;
constexpr double protonMassC2 = protonMass * cSquared;
constexpr double neutronMass = 1.67492749804e-24 * gram;
constexpr double neutronMassC2 = neutronMass * cSquared;
constexpr double electronMass = 9.1093837015e-28 * gram;
constexpr double electronMassC2 = electronMass * cSquared;
constexpr double sunMass = 1.989e33 * gram;
constexpr double hPlanck = 6.62607015e-34 * joule * second;   // exact
constexpr double kBoltzmann = 1.380649e-23 * joule / kelvin;  // CODATA 2018 exact
constexpr double electronRadius = 2.8179403262e-15 * meter;
constexpr double IsH = 19 * eV;   // H  eff. ioniz. potential
constexpr double IsHe = 44 * eV;  // He eff. ioniz. potential
constexpr double barn = 1e-24 * cm2;
constexpr double mbarn = 1e-3 * barn;

// MOMENTUM UNITS
constexpr double eV_c = electronvolt / cLight;
constexpr double keV_c = 1e3 * eV_c;
constexpr double MeV_c = 1e6 * eV_c;
constexpr double GeV_c = 1e9 * eV_c;
constexpr double TeV_c = 1e12 * eV_c;
constexpr double PeV_c = 1e15 * eV_c;

// CODE-SPECIFIC CONSTANTS
constexpr double E_SN = 1e51 * erg;
constexpr double snRate = 1. / (30. * year);
constexpr double galaxySize = 10. * kpc;
constexpr double f_He = 0.08;
constexpr double K_He = 2.51984209979;  // 4^(2/3)
constexpr double meanISMmass = protonMass * (1. + 4 * f_He) / (1. + f_He);
constexpr double inelasticity = 0.5;

}  // namespace CGS
}  // namespace CRAMS

#endif  // CRAMS_CORE_CGS_H_
