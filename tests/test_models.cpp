// Verifies that every shipped cross-section model file exists and is formatted
// per the hard-coded parameters. Constructing a model reads its CSV and
// validates the first-row energy grid against the grid built from
// T_min/T_max/T_size and checks the per-row column counts, throwing on any
// mismatch — so a successful construction *is* the format check.
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

#include "crams/core/cgs.h"
#include "crams/core/pid.h"
#include "crams/fragmentation.h"
#include "crams/inelastic.h"

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

// Returns true if *construct* runs without throwing (file present + well formatted).
static bool loads(const std::function<void()>& construct) {
  try {
    construct();
  } catch (const std::exception& e) {
    std::cerr << "  model load failed: " << e.what() << "\n";
    return false;
  }
  return true;
}

void test_inelastic_models_load() {
  CHECK(loads([] { CRAMS::InXsecTripathi99 m; }));
  CHECK(loads([] { CRAMS::InXsecGlauber m; }));
  CHECK(loads([] { CRAMS::InXsecCrosec m; }));
}

void test_fragmentation_models_load() {
  CHECK(loads([] { CRAMS::NucFragFluka4Dragon m; }));
  CHECK(loads([] { CRAMS::NucFragUsineGalprop17Opt12 m; }));
  CHECK(loads([] { CRAMS::NucFragUsineGalprop17Opt22 m; }));
  CHECK(loads([] { CRAMS::NucFragUsineWebber03Coste12 m; }));
}

void test_inelastic_values_sane() {
  const double T = 10. * CRAMS::CGS::GeV;
  CRAMS::InXsecTripathi99 trip;
  CRAMS::InXsecGlauber glau;
  CRAMS::InXsecCrosec cro;
  const CRAMS::InelasticXsec* models[] = {&trip, &glau, &cro};
  for (const CRAMS::InelasticXsec* m : models) {
    CHECK(m->getXsecOnHtarget(CRAMS::C12, T) > 0.);                                   // table populated
    CHECK(m->getXsecOnHtarget(CRAMS::Fe56, T) > m->getXsecOnHtarget(CRAMS::C12, T));  // heavier -> larger
  }
}

void test_fragmentation_values_sane() {
  const double T = 10. * CRAMS::CGS::GeV;
  CRAMS::NucFragFluka4Dragon fluka;
  // C12 -> B11 is a major fragmentation channel; must be positive and finite.
  const double sigma = fluka.getXsecOnHtarget(CRAMS::C12, CRAMS::B11, T);
  CHECK(sigma > 0.);
}

// Minimal table model with the inelastic grid (448 pts), pointed at an arbitrary
// file, so we can feed it a malformed table and confirm the loader rejects it.
namespace {
class TestInelasticTable : public CRAMS::InXsecFromTable {
 public:
  explicit TestInelasticTable(const std::string& file)
      : CRAMS::InXsecFromTable("TestModel", file, 0.01 * CRAMS::CGS::GeV, 1e5 * CRAMS::CGS::GeV, 448) {}
};
}  // namespace

void test_malformed_grid_is_rejected() {
  // Header advertises only 2 grid columns, but the model's hard-coded grid is
  // 448 points -> construction must throw.
  const std::string path = "/tmp/test_models_badgrid.csv";
  std::ofstream(path) << "Z,A,0.01,0.02\n1,2,100.0,100.0\n";
  bool threw = false;
  try {
    TestInelasticTable m(path);
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);
}

int main() {
  test_inelastic_models_load();
  test_fragmentation_models_load();
  test_inelastic_values_sane();
  test_fragmentation_values_sane();
  test_malformed_grid_is_rejected();

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
