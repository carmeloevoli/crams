#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

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

#define CHECK_THROW(expr, exc) \
  do {                         \
    bool caught_ = false;      \
    try {                      \
      (void)(expr);            \
    } catch (const exc&) {     \
      caught_ = true;          \
    }                          \
    CHECK(caught_);            \
  } while (0)

static bool approx(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol * std::abs(b) + tol; }

static void write_file(const std::string& path, const std::string& content) {
  std::ofstream f(path);
  assert(f.is_open());
  f << content;
}

static size_t count_nucleilist_rows() {
  std::ifstream f("data/nucleilist.csv");
  assert(f.is_open());

  size_t count = 0;
  std::string line;
  while (std::getline(f, line)) {
    if (!line.empty() && line[0] != '#') ++count;
  }
  return count;
}

void test_default_list_and_injection_parameters() {
  CRAMS::ParticleList particles;
  const auto& list = particles.getList();

  CHECK(list.size() == count_nucleilist_rows());
  CHECK(list.find(CRAMS::H1) != list.end());
  CHECK(list.find(CRAMS::He4) != list.end());
  CHECK(list.find(CRAMS::C12) != list.end());
  CHECK(list.find(CRAMS::C14) != list.end());

  CHECK(approx(list.at(CRAMS::H1).abundance, 5.06605e-02 * 0.9999806));
  CHECK(approx(list.at(CRAMS::H2).abundance, 5.06605e-02 * 0.0000194));
  CHECK(approx(list.at(CRAMS::He4).abundance, 2.54369e-02 * 0.99983403));
  CHECK(approx(list.at(CRAMS::C12).abundance, 3.98879e-03 * 0.988922));

  CHECK(approx(list.at(CRAMS::H1).slope, 4.37486));
  CHECK(approx(list.at(CRAMS::He4).slope, 4.30995));
  CHECK(approx(list.at(CRAMS::C12).slope, 4.32798));
}

void test_read_params_updates_charge_groups() {
  const std::string path = "/tmp/test_particlelist_params.ini";
  write_file(path, "# ignored\nbadline\nqC 1.0e-2\nH_slope 4.1\nHe_slope 4.2\nslope 4.9\n");

  CRAMS::ParticleList particles;
  particles.readParamsFromFile(path);
  const auto& list = particles.getList();

  CHECK(approx(list.at(CRAMS::C12).abundance, 1.0e-2 * 0.988922));
  CHECK(approx(list.at(CRAMS::C13).abundance, 1.0e-2 * 0.011078));
  CHECK(approx(list.at(CRAMS::C14).abundance, 0.0));

  CHECK(approx(list.at(CRAMS::H1).slope, 4.1));
  CHECK(approx(list.at(CRAMS::He4).slope, 4.2));
  CHECK(approx(list.at(CRAMS::C12).slope, 4.9));
}

void test_read_params_missing_file_throws() {
  CRAMS::ParticleList particles;
  CHECK_THROW(particles.readParamsFromFile("/tmp/no_such_particlelist_file.ini"), std::runtime_error);
}

void test_copy_rebuilds_charge_index() {
  const std::string path = "/tmp/test_particlelist_copy_params.ini";
  write_file(path, "qC 2.0e-2\n");

  CRAMS::ParticleList original;
  CRAMS::ParticleList copy = original;
  copy.readParamsFromFile(path);

  CHECK(approx(copy.getList().at(CRAMS::C12).abundance, 2.0e-2 * 0.988922));
  CHECK(approx(original.getList().at(CRAMS::C12).abundance, 3.98879e-03 * 0.988922));
}

void test_move_rebuilds_charge_index() {
  const std::string path = "/tmp/test_particlelist_move_params.ini";
  write_file(path, "qHe 3.0e-2\n");

  CRAMS::ParticleList original;
  CRAMS::ParticleList moved = std::move(original);
  moved.readParamsFromFile(path);

  CHECK(approx(moved.getList().at(CRAMS::He3).abundance, 3.0e-2 * 0.00016597));
  CHECK(approx(moved.getList().at(CRAMS::He4).abundance, 3.0e-2 * 0.99983403));
}

void test_mutable_getList_keeps_charge_updates_safe() {
  const std::string path = "/tmp/test_particlelist_mutable_params.ini";
  write_file(path, "qC 4.0e-2\n");

  CRAMS::ParticleList particles;
  auto& list = particles.getList();
  const CRAMS::PID C15{6, 15};
  list.erase(CRAMS::C13);
  list.emplace(C15, CRAMS::NucleusParameters{0., 4.0, 0.5, -1., true, false});

  particles.readParamsFromFile(path);

  CHECK(list.find(CRAMS::C13) == list.end());
  CHECK(approx(list.at(CRAMS::C12).abundance, 4.0e-2 * 0.988922));
  CHECK(approx(list.at(C15).abundance, 4.0e-2 * 0.5));
}

int main() {
  test_default_list_and_injection_parameters();
  test_read_params_updates_charge_groups();
  test_read_params_missing_file_throws();
  test_copy_rebuilds_charge_index();
  test_move_rebuilds_charge_index();
  test_mutable_getList_keeps_charge_updates_safe();

  if (g_fail != 0) {
    std::cerr << g_fail << " particlelist test(s) failed, " << g_pass << " passed\n";
    return 1;
  }

  std::cout << g_pass << " particlelist tests passed\n";
  return 0;
}
