#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>

#include "crams/core/input.h"
#include "crams/fragmentation.h"
#include "crams/inelastic.h"
#include "crams/particle.h"
#include "crams/particlelist.h"

static int g_pass = 0;
static int g_fail = 0;

#define CHECK(cond)                                                                    \
  do {                                                                                 \
    if (cond) {                                                                        \
      ++g_pass;                                                                        \
    } else {                                                                           \
      ++g_fail;                                                                        \
      std::cerr << "FAIL: " << #cond << " at " << __FILE__ << ":" << __LINE__ << "\n"; \
    }                                                                                  \
  } while (0)

static bool approx(double a, double b, double tol = 1e-12) { return std::abs(a - b) <= tol * std::abs(b) + tol; }

void write_file(const std::string& path, const std::string& contents) {
  std::ofstream out(path);
  out << contents;
}

void test_particle_ownership_traits() {
  CHECK(std::is_move_constructible<CRAMS::Particle>::value);
  CHECK(std::is_move_assignable<CRAMS::Particle>::value);
  CHECK(!std::is_copy_constructible<CRAMS::Particle>::value);
  CHECK(!std::is_copy_assignable<CRAMS::Particle>::value);
}

void test_constructor_copies_nucleus_parameters() {
  const CRAMS::NucleusParameters params{1.5, 4.2, 0.3, 7.0, false, true};
  CRAMS::Particle particle(CRAMS::C12, params);

  CHECK(particle.getPid() == CRAMS::C12);
  CHECK(approx(particle.getAbundance(), 1.5));
  CHECK(approx(particle.getSlope(), 4.2));
  CHECK(approx(particle.getDecayTime(), 7.0));
  CHECK(!particle.isStable());
}

void test_build_vectors_initializes_intensity() {
  CRAMS::Input input;
  CRAMS::Particle particle(CRAMS::H1);
  particle.buildVectors(input);

  CHECK(particle.getEnergyVector().size() == input.TSimSize());
  CHECK(particle.getIntensityVector().size() == input.TSimSize());
  CHECK(particle.getIntensityVector().front() == 0.);
  CHECK(particle.getIntensityVector().back() == 0.);
}

void test_empty_interpolation_is_zero() {
  CRAMS::Particle particle(CRAMS::H1);
  CHECK(particle.I_T_interpol(1.0) == 0.);
}

void test_move_preserves_particle_identity() {
  CRAMS::Particle particle(CRAMS::He4);
  CRAMS::Particle moved(std::move(particle));
  CHECK(moved.getPid() == CRAMS::He4);
}

void test_crank_nicolson_solver_computes_primary_flux() {
  const std::string path = "/tmp/test_particle_crank_nicolson.ini";
  write_file(path, "solver crank_nicolson\n");

  CRAMS::Input input;
  input.readParamsFromFile(path);

  const CRAMS::NucleusParameters protonParams{1., 4.2, 1., -1., true, true};
  CRAMS::Particle particle(CRAMS::H1, protonParams);
  CRAMS::InXsecTripathi99 inelasticXsecs;
  CRAMS::NucFragFluka4Dragon nucfragXsecs;
  const CRAMS::Particles noParents;

  particle.buildVectors(input);
  particle.buildGrammage(input);
  particle.buildLosses(input);
  particle.buildPrimarySource(input);
  particle.buildInelasticXsecs(inelasticXsecs);
  particle.buildSecondarySource(input, noParents, nucfragXsecs);
  particle.computeIntensity(input);

  const auto& intensity = particle.getIntensityVector();
  CHECK(particle.isDone());
  CHECK(intensity.size() == input.TSimSize());
  CHECK(intensity.front() > 0.);
  CHECK(intensity.back() == 0.);
  CHECK(std::all_of(intensity.begin(), intensity.end(),
                    [](double value) { return std::isfinite(value) && value >= 0.; }));
}

void test_compute_intensity_without_secondary_or_inelastic_sources() {
  const std::string path = "/tmp/test_particle_no_optional_sources.ini";
  write_file(path, "solver crank_nicolson\n");

  CRAMS::Input input;
  input.readParamsFromFile(path);

  const CRAMS::NucleusParameters protonParams{1., 4.2, 1., -1., true, true};
  CRAMS::Particle particle(CRAMS::H1, protonParams);
  CRAMS::NucFragFluka4Dragon nucfragXsecs;
  const CRAMS::Particles noParents;

  particle.buildVectors(input);
  particle.buildGrammage(input);
  particle.buildLosses(input);
  particle.buildPrimarySource(input);
  particle.buildSecondarySource(input, noParents, nucfragXsecs);
  particle.reset();

  particle.buildVectors(input);
  particle.buildGrammage(input);
  particle.buildLosses(input);
  particle.buildPrimarySource(input);
  particle.computeIntensity(input);

  const auto& intensity = particle.getIntensityVector();
  CHECK(particle.isDone());
  CHECK(intensity.size() == input.TSimSize());
  CHECK(intensity.front() > 0.);
  CHECK(intensity.back() == 0.);
  CHECK(std::all_of(intensity.begin(), intensity.end(),
                    [](double value) { return std::isfinite(value) && value >= 0.; }));
}

void test_dump_without_secondary_or_inelastic_sources() {
  CRAMS::Input input;
  const CRAMS::NucleusParameters protonParams{1., 4.2, 1., -1., true, true};
  CRAMS::Particle particle(CRAMS::H1, protonParams);

  particle.buildVectors(input);
  particle.buildGrammage(input);
  particle.buildLosses(input);
  particle.buildPrimarySource(input);

  bool dumped = false;
  try {
    particle.dump();
    dumped = true;
  } catch (const std::exception&) {
  }

  std::ifstream dumpFile("output/crams_particle_dump_1_1.txt");
  std::string contents;
  std::string firstRow;
  for (std::string line; std::getline(dumpFile, line);) {
    contents += line + "\n";
    if (!line.empty() && line[0] != '#') {
      firstRow = line;
      break;
    }
  }

  CHECK(dumped);
  CHECK(contents.find("# T [GeV] -> 1") != std::string::npos);
  CHECK(contents.find("# Q_sec -> 4") != std::string::npos);
  CHECK(contents.find("# tau_inelastic [Myr] -> 11") != std::string::npos);
  CHECK(!firstRow.empty());
}

int main() {
  test_particle_ownership_traits();
  test_constructor_copies_nucleus_parameters();
  test_build_vectors_initializes_intensity();
  test_empty_interpolation_is_zero();
  test_move_preserves_particle_identity();
  test_crank_nicolson_solver_computes_primary_flux();
  test_compute_intensity_without_secondary_or_inelastic_sources();
  test_dump_without_secondary_or_inelastic_sources();

  if (g_fail != 0) {
    std::cerr << g_fail << " particle test(s) failed, " << g_pass << " passed\n";
    return 1;
  }

  std::cout << g_pass << " particle tests passed\n";
  return 0;
}
