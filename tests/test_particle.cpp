#include <cmath>
#include <iostream>
#include <type_traits>

#include "crams/particle.h"
#include "crams/particlelist.h"
#include "crams/core/input.h"

static int g_pass = 0;
static int g_fail = 0;

#define CHECK(cond)                                                                      \
  do {                                                                                   \
    if (cond) {                                                                          \
      ++g_pass;                                                                          \
    } else {                                                                             \
      ++g_fail;                                                                          \
      std::cerr << "FAIL: " << #cond << " at " << __FILE__ << ":" << __LINE__ << "\n"; \
    }                                                                                    \
  } while (0)

static bool approx(double a, double b, double tol = 1e-12) {
  return std::abs(a - b) <= tol * std::abs(b) + tol;
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

int main() {
  test_particle_ownership_traits();
  test_constructor_copies_nucleus_parameters();
  test_build_vectors_initializes_intensity();
  test_empty_interpolation_is_zero();
  test_move_preserves_particle_identity();

  if (g_fail != 0) {
    std::cerr << g_fail << " particle test(s) failed, " << g_pass << " passed\n";
    return 1;
  }

  std::cout << g_pass << " particle tests passed\n";
  return 0;
}
