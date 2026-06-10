#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "crams/core/input.h"
#include "crams/core/output.h"
#include "crams/particle.h"

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

std::string read_file(const std::string& path) {
  std::ifstream in(path);
  std::ostringstream contents;
  contents << in.rdbuf();
  return contents.str();
}

void test_spectra_headers_map_species_to_columns() {
  CRAMS::Input input;
  CRAMS::Particles particles;
  particles.emplace_back(CRAMS::H1);
  particles.emplace_back(CRAMS::He3);

  CRAMS::OutputManager output(particles, input);
  output.dumpSpectraRigidity();
  output.dumpSpectraEkn();
  output.dumpIsotopes();

  const auto rigidity = read_file("output/test_spectra_R_0.txt");
  const auto ekn = read_file("output/test_spectra_Ekn_0.txt");
  const auto isotopes = read_file("output/test_isotopes_R_0.txt");

  CHECK(rigidity.find("# R [GV] -> 1") != std::string::npos);
  CHECK(rigidity.find("# H -> 2") != std::string::npos);
  CHECK(rigidity.find("# He -> 3") != std::string::npos);
  CHECK(rigidity.find("# H1 -> 2") == std::string::npos);
  CHECK(rigidity.find("# Flux unit -> 1 / (GeV m2 s sr)") != std::string::npos);
  CHECK(rigidity.find("pbar") == std::string::npos);

  CHECK(ekn.find("# T [GeV/n] -> 1") != std::string::npos);
  CHECK(ekn.find("# H -> 2") != std::string::npos);
  CHECK(ekn.find("# He -> 3") != std::string::npos);

  CHECK(isotopes.find("# Be9 -> 2") != std::string::npos);
  CHECK(isotopes.find("# Be10 -> 3") != std::string::npos);
}

int main() {
  test_spectra_headers_map_species_to_columns();

  if (g_fail != 0) {
    std::cerr << g_fail << " output test(s) failed, " << g_pass << " passed\n";
    return 1;
  }

  std::cout << g_pass << " output tests passed\n";
  return 0;
}
