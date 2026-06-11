#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "crams/utils/grid.h"

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

void test_default_constructor() {
  CRAMS::Grid<double> g;
  CHECK(g.getNx() == 0);
  CHECK(g.getNy() == 0);
  CHECK(g.size() == 0);
  CHECK(g.get().empty());
}

void test_sized_constructor() {
  CRAMS::Grid<double> g(4, 3);
  CHECK(g.getNx() == 4);
  CHECK(g.getNy() == 3);
  CHECK(g.size() == 12);
  CHECK(g.get().size() == 12);
}

void test_index_layout() {
  // index(ix, iy) = iy + Ny * ix — y is the fast index (column-major in y)
  CRAMS::Grid<double> g(3, 4);  // Nx=3, Ny=4
  CHECK(g.index(0, 0) == 0);
  CHECK(g.index(0, 1) == 1);
  CHECK(g.index(0, 3) == 3);
  CHECK(g.index(1, 0) == 4);   // iy=0 + Ny*ix=4*1
  CHECK(g.index(2, 3) == 11);  // iy=3 + 4*2 = 11
}

void test_get_set_2d() {
  CRAMS::Grid<double> g(3, 4);
  g.get(0, 0) = 1.0;
  g.get(1, 2) = 2.5;
  g.get(2, 3) = -7.0;
  CHECK(g.get(0, 0) == 1.0);
  CHECK(g.get(1, 2) == 2.5);
  CHECK(g.get(2, 3) == -7.0);
}

void test_get_set_const() {
  CRAMS::Grid<double> g(2, 2);
  g.get(0, 1) = 42.0;
  const auto& cg = g;
  CHECK(cg.get(0, 1) == 42.0);
  CHECK(cg.get().size() == 4);
}

void test_flat_vector_accessor() {
  CRAMS::Grid<int> g(2, 3);
  auto& v = g.get();
  for (int i = 0; i < 6; ++i) v[i] = i;
  CHECK(g.get(0, 0) == 0);
  CHECK(g.get(0, 1) == 1);
  CHECK(g.get(1, 0) == 3);
  CHECK(g.get(1, 2) == 5);
}

void test_copy_lvalue() {
  CRAMS::Grid<double> g(2, 3);
  std::vector<double> data = {1., 2., 3., 4., 5., 6.};
  g.copy(data);
  CHECK(g.get()[0] == 1.0);
  CHECK(g.get()[5] == 6.0);
  // original data unmodified (copy, not move)
  CHECK(data.size() == 6);
}

void test_copy_rvalue() {
  CRAMS::Grid<double> g(2, 3);
  std::vector<double> data = {10., 20., 30., 40., 50., 60.};
  g.copy(std::move(data));
  CHECK(g.get()[0] == 10.0);
  CHECK(g.get()[5] == 60.0);
}

void test_copy_preserves_2d_indexing() {
  // copy loads flat data; verify 2d indexing matches
  // index(ix, iy) = iy + Ny*ix, so for Nx=2, Ny=3:
  // (0,0)=0, (0,1)=1, (0,2)=2, (1,0)=3, (1,1)=4, (1,2)=5
  CRAMS::Grid<double> g(2, 3);
  g.copy({0., 1., 2., 3., 4., 5.});
  CHECK(g.get(0, 0) == 0.);
  CHECK(g.get(0, 2) == 2.);
  CHECK(g.get(1, 0) == 3.);
  CHECK(g.get(1, 2) == 5.);
}

void test_for_each_scale() {
  CRAMS::Grid<double> g(2, 3);
  g.copy({1., 2., 3., 4., 5., 6.});
  g.for_each([](double& v) { v *= 2.; });
  CHECK(g.get()[0] == 2.0);
  CHECK(g.get()[5] == 12.0);
}

void test_for_each_lambda_capture() {
  CRAMS::Grid<double> g(3, 3);
  g.copy({1., 2., 3., 4., 5., 6., 7., 8., 9.});
  const double factor = 3.0;
  g.for_each([factor](double& v) { v *= factor; });
  CHECK(g.get()[0] == 3.0);
  CHECK(g.get()[8] == 27.0);
}

void test_for_each_count_elements() {
  CRAMS::Grid<int> g(4, 5);
  int count = 0;
  g.for_each([&count](int&) { ++count; });
  CHECK(count == 20);
}

void test_min_max_basic() {
  CRAMS::Grid<double> g(2, 3);
  g.copy({3., 1., 4., 1., 5., 9.});
  CHECK(g.min() == 1.0);
  CHECK(g.max() == 9.0);
}

void test_min_max_uniform() {
  CRAMS::Grid<double> g(3, 3);
  g.copy({7., 7., 7., 7., 7., 7., 7., 7., 7.});
  CHECK(g.min() == 7.0);
  CHECK(g.max() == 7.0);
}

void test_min_max_negative() {
  CRAMS::Grid<double> g(1, 4);
  g.copy({-5., -2., 0., 3.});
  CHECK(g.min() == -5.0);
  CHECK(g.max() == 3.0);
}

void test_int_grid() {
  CRAMS::Grid<int> g(3, 2);
  g.copy({10, 20, 30, 40, 50, 60});
  CHECK(g.min() == 10);
  CHECK(g.max() == 60);
  CHECK(g.get(1, 0) == 30);
}

void test_single_element_grid() {
  CRAMS::Grid<double> g(1, 1);
  CHECK(g.size() == 1);
  g.get(0, 0) = 99.0;
  CHECK(g.min() == 99.0);
  CHECK(g.max() == 99.0);
}

void test_large_grid_iota() {
  const size_t nx = 100, ny = 50;
  CRAMS::Grid<double> g(nx, ny);
  auto& v = g.get();
  std::iota(v.begin(), v.end(), 0.0);
  CHECK(g.size() == nx * ny);
  CHECK(g.min() == 0.0);
  CHECK(g.max() == static_cast<double>(nx * ny - 1));
}

void test_for_each_sum() {
  CRAMS::Grid<double> g(3, 4);
  g.copy({1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.});
  double sum = 0.;
  g.for_each([&sum](double& v) { sum += v; });
  CHECK(std::abs(sum - 12.) < 1e-12);
}

int main() {
  test_default_constructor();
  test_sized_constructor();
  test_index_layout();
  test_get_set_2d();
  test_get_set_const();
  test_flat_vector_accessor();
  test_copy_lvalue();
  test_copy_rvalue();
  test_copy_preserves_2d_indexing();
  test_for_each_scale();
  test_for_each_lambda_capture();
  test_for_each_count_elements();
  test_min_max_basic();
  test_min_max_uniform();
  test_min_max_negative();
  test_int_grid();
  test_single_element_grid();
  test_large_grid_iota();
  test_for_each_sum();

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
