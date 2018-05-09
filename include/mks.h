#ifndef UNITS_H_
#define UNITS_H_

#include <cmath>

// powers
#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define pow4(A) ((A)*(A)*(A)*(A))
#define pow5(A) ((A)*(A)*(A)*(A)*(A))

// MKS units
static const double meter = 1;
static const double kilogram = 1;
static const double second = 1;
static const double ampere = 1;
static const double kelvin = 1;
static const double mole = 1;
static const double joule = kilogram * pow2(meter) / pow2(second);
static const double tesla = kilogram / ampere / pow2(second);
static const double coulomb = ampere * second;
static const double sr = 1;

// derived units
static const double centimeter = 1e-2 * meter;
static const double kilometer = 1e3 * meter;
static const double gram = 1e-3 * kilogram;
static const double mgram = 1e-6 * kilogram;
static const double erg = 1e-7 * joule;
static const double gauss = 1e-4 * tesla;
static const double microgauss = 1e-6 * gauss;
static const double nanogauss = 1e-9 * gauss;
static const double barn = 1e-28 * pow2(meter);
static const double mbarn = 1e-3 * barn;

// abbreviations
static const double sec = second;
static const double K = kelvin;
static const double kg = kilogram;
static const double cm = centimeter;
static const double km = kilometer;
static const double cm2 = cm * cm;
static const double cm3 = cm * cm * cm;
static const double m2 = meter * meter;
static const double m3 = meter * meter * meter;
static const double muG = microgauss;
static const double nG = nanogauss;

// electron volt
static const double electronvolt = 1.60217657e-19 * joule;
static const double kiloelectronvolt = 1e3 * electronvolt;
static const double megaelectronvolt = 1e6 * electronvolt;
static const double gigaelectronvolt = 1e9 * electronvolt;
static const double teraelectronvolt = 1e12 * electronvolt;
static const double petaelectronvolt = 1e15 * electronvolt;
static const double exaelectronvolt = 1e18 * electronvolt;
static const double eV = electronvolt;
static const double keV = kiloelectronvolt;
static const double MeV = megaelectronvolt;
static const double GeV = gigaelectronvolt;
static const double TeV = teraelectronvolt;
static const double PeV = petaelectronvolt;
static const double EeV = exaelectronvolt;

// time
static const double year = 3.15569e7 * second;
static const double kiloyear = 1e3 * year;
static const double megayear = 1e6 * year;
static const double gigayear = 1e9 * year;
static const double kyr = kiloyear;
static const double Myr = megayear;
static const double Gyr = gigayear;

// parsec
static const double parsec = 3.0856775807e18 * centimeter;
static const double kiloparsec = 1e3 * parsec;
static const double megaparsec = 1e6 * parsec;
static const double gigaparsec = 1e9 * parsec;
static const double pc = parsec;
static const double kpc = kiloparsec;
static const double Mpc = megaparsec;
static const double Gpc = gigaparsec;

// physical constants
static const double c_light = 2.99792458e10 * centimeter / second;
static const double c_squared = c_light * c_light;
static const double electron_mass = 9.10938356e-31 * kilogram;
static const double electron_mass_c = electron_mass * c_light;
static const double electron_mass_c2 = electron_mass * c_squared;
static const double proton_mass = 1.6726219e-27 * kilogram;
static const double proton_mass_c = proton_mass * c_light;
static const double proton_mass_c2 = proton_mass * c_squared;
static const double sun_mass = 1.989e30 * kg;
static const double sigma_th = 6.6524e-29 * pow2(meter);
static const double h_planck = 6.62607004e-34 * pow2(meter) * kg / sec;
static const double k_boltzmann = 1.3806488e-23 * joule / kelvin;
static const double elementary_charge = 1.60217662e-19 * coulomb;
static const double vacuum_permeability = 4e-7 * M_PI * tesla * meter / ampere;
static const double electron_radius = 2.8179403227e-15 * meter;

// derived units with physical constants
static const double keV_c = keV / c_light;
static const double MeV_c = MeV / c_light;
static const double GeV_c = GeV / c_light;
static const double TeV_c = TeV / c_light;
static const double PeV_c = PeV / c_light;
static const double EeV_c = EeV / c_light;

#endif
