#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "crams/core/cgs.h"
#include "crams/core/input.h"

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

// --- default values ---

void test_default_TSimMin() { CHECK(approx(CRAMS::Input{}.TSimMin(), 0.1 * CRAMS::CGS::GeV)); }
void test_default_TSimMax() { CHECK(approx(CRAMS::Input{}.TSimMax(), 100. * CRAMS::CGS::TeV)); }
void test_default_TSimSize() { CHECK(CRAMS::Input{}.TSimSize() == 300); }
void test_default_ROutputMin() { CHECK(approx(CRAMS::Input{}.ROutputMin(), 1. * CRAMS::CGS::GeV)); }
void test_default_ROutputMax() { CHECK(approx(CRAMS::Input{}.ROutputMax(), 10. * CRAMS::CGS::TeV)); }
void test_default_ROutputSize() { CHECK(CRAMS::Input{}.ROutputSize() == 100); }
void test_default_doSecondary() { CHECK(CRAMS::Input{}.doSecondary() == true); }
void test_default_H() { CHECK(approx(CRAMS::Input{}.H(), 7. * CRAMS::CGS::kpc)); }
void test_default_mu() { CHECK(approx(CRAMS::Input{}.mu(), 2.3 * CRAMS::CGS::mgram / CRAMS::CGS::cm2)); }
void test_default_delta() { CHECK(approx(CRAMS::Input{}.delta(), 5.65132e-01)); }
void test_default_ddelta() { CHECK(approx(CRAMS::Input{}.ddelta(), 0.22)); }
void test_default_phi() { CHECK(approx(CRAMS::Input{}.modulationPotential(), 4.87754e-01 * CRAMS::CGS::GeV)); }
void test_default_X_s_negative() { CHECK(CRAMS::Input{}.X_s() < 0.); }
void test_default_simname() { CHECK(CRAMS::Input{}.simname() == "test"); }
void test_default_id() { CHECK(CRAMS::Input{}.id() == 0); }
void test_default_flux_solver() {
  CHECK(CRAMS::Input{}.fluxSolver() == CRAMS::FluxSolver::CrankNicolson);
  CHECK(CRAMS::Input{}.fluxSolverName() == "crank_nicolson");
}

// --- setParam ---

void test_setParam_D0() {
  CRAMS::Input in;
  in.setParam("D0", 3.0);
  CHECK(approx(in.D_0(), 3.0 * 1e28 * CRAMS::CGS::cm2 / CRAMS::CGS::sec));
}

void test_setParam_H() {
  CRAMS::Input in;
  in.setParam("H", 5.0);
  CHECK(approx(in.H(), 5.0 * CRAMS::CGS::kpc));
}

void test_setParam_delta() {
  CRAMS::Input in;
  in.setParam("delta", 0.45);
  CHECK(approx(in.delta(), 0.45));
}

void test_setParam_ddelta() {
  CRAMS::Input in;
  in.setParam("ddelta", 0.1);
  CHECK(approx(in.ddelta(), 0.1));
}

void test_setParam_Rb() {
  CRAMS::Input in;
  in.setParam("Rb", 300.0);
  CHECK(approx(in.R_b(), 300.0 * CRAMS::CGS::GeV));
}

void test_setParam_vA() {
  CRAMS::Input in;
  in.setParam("vA", 10.0);
  CHECK(approx(in.v_A(), 10.0 * CRAMS::CGS::km / CRAMS::CGS::sec));
}

void test_setParam_phi() {
  CRAMS::Input in;
  in.setParam("phi", 0.6);
  CHECK(approx(in.modulationPotential(), 0.6 * CRAMS::CGS::GeV));
}

void test_setParam_xs() {
  CRAMS::Input in;
  in.setParam("xs", 0.5);
  CHECK(approx(in.X_s(), 0.5 * CRAMS::CGS::gram / CRAMS::CGS::cm2));
}

void test_setParam_id() {
  CRAMS::Input in;
  in.setParam("id", 7.0);
  CHECK(in.id() == 7);
}

void test_setParam_key_case_insensitive() {
  CRAMS::Input in;
  in.setParam("DELTA", 0.3);
  CHECK(approx(in.delta(), 0.3));
}

void test_setParam_unknown_key_ignored() {
  CRAMS::Input in;
  const double prev = in.delta();
  in.setParam("unknownparam", 99.0);
  CHECK(approx(in.delta(), prev));
}

// --- setSimname ---

void test_setSimname_strips_ini() {
  CRAMS::Input in;
  in.setSimname("params.ini");
  CHECK(in.simname() == "params");
}

void test_setSimname_path_uses_basename() {
  CRAMS::Input in;
  in.setSimname("run/my_sim.ini");
  CHECK(in.simname() == "my_sim");
}

void test_setSimname_wrong_extension_throws() {
  CRAMS::Input in;
  CHECK_THROW(in.setSimname("params.txt"), std::runtime_error);
}

void test_setSimname_no_extension_throws() {
  CRAMS::Input in;
  CHECK_THROW(in.setSimname("params"), std::runtime_error);
}

// --- readParamsFromFile ---

void test_readParamsFromFile_basic() {
  const std::string path = "/tmp/test_input_basic.ini";
  write_file(path, "H 10.0\ndelta 0.4\n");
  CRAMS::Input in;
  in.readParamsFromFile(path);
  CHECK(approx(in.H(), 10.0 * CRAMS::CGS::kpc));
  CHECK(approx(in.delta(), 0.4));
}

void test_readParamsFromFile_ignores_bad_lines() {
  const std::string path = "/tmp/test_input_badlines.ini";
  write_file(path, "# comment\ndelta 0.5\n\nbadline\nH 3.0\n");
  CRAMS::Input in;
  in.readParamsFromFile(path);
  CHECK(approx(in.delta(), 0.5));
  CHECK(approx(in.H(), 3.0 * CRAMS::CGS::kpc));
}

void test_readParamsFromFile_missing_file_throws() {
  CRAMS::Input in;
  CHECK_THROW(in.readParamsFromFile("/tmp/no_such_file_xyzzy.ini"), std::runtime_error);
}

void test_readParamsFromFile_solver_string() {
  const std::string path = "/tmp/test_input_solver.ini";
  write_file(path, "solver exponential\n");
  CRAMS::Input in;
  in.readParamsFromFile(path);
  CHECK(in.fluxSolver() == CRAMS::FluxSolver::Exponential);
}

void test_readParamsFromFile_solver_string_with_underscore() {
  const std::string path = "/tmp/test_input_solver_underscore.ini";
  write_file(path, "solver crank_nicolson\n");
  CRAMS::Input in;
  in.readParamsFromFile(path);
  CHECK(in.fluxSolver() == CRAMS::FluxSolver::CrankNicolson);
}

void test_readParamsFromFile_solver_invalid_throws() {
  const std::string path = "/tmp/test_input_solver_invalid.ini";
  write_file(path, "solver mystery\n");
  CRAMS::Input in;
  CHECK_THROW(in.readParamsFromFile(path), std::runtime_error);
}

// --- copyability ---

void test_input_is_copyable() {
  CRAMS::Input a;
  a.setParam("H", 4.0);
  CRAMS::Input b = a;
  CHECK(approx(b.H(), 4.0 * CRAMS::CGS::kpc));
  // modifying copy does not affect original
  b.setParam("H", 8.0);
  CHECK(approx(a.H(), 4.0 * CRAMS::CGS::kpc));
  CHECK(approx(b.H(), 8.0 * CRAMS::CGS::kpc));
}

int main() {
  test_default_TSimMin();
  test_default_TSimMax();
  test_default_TSimSize();
  test_default_ROutputMin();
  test_default_ROutputMax();
  test_default_ROutputSize();
  test_default_doSecondary();
  test_default_H();
  test_default_mu();
  test_default_delta();
  test_default_ddelta();
  test_default_phi();
  test_default_X_s_negative();
  test_default_simname();
  test_default_id();
  test_default_flux_solver();

  test_setParam_D0();
  test_setParam_H();
  test_setParam_delta();
  test_setParam_ddelta();
  test_setParam_Rb();
  test_setParam_vA();
  test_setParam_phi();
  test_setParam_xs();
  test_setParam_id();
  test_setParam_key_case_insensitive();
  test_setParam_unknown_key_ignored();

  test_setSimname_strips_ini();
  test_setSimname_path_uses_basename();
  test_setSimname_wrong_extension_throws();
  test_setSimname_no_extension_throws();

  test_readParamsFromFile_basic();
  test_readParamsFromFile_ignores_bad_lines();
  test_readParamsFromFile_missing_file_throws();
  test_readParamsFromFile_solver_string();
  test_readParamsFromFile_solver_string_with_underscore();
  test_readParamsFromFile_solver_invalid_throws();

  test_input_is_copyable();

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
