#ifndef CRAMS_UTILS_GRID_H_
#define CRAMS_UTILS_GRID_H_

#include <algorithm>
#include <cassert>
#include <vector>

namespace CRAMS {

template <typename T>
class Grid {
 public:
  Grid() = default;
  Grid(size_t nx, size_t ny) : m_Nx(nx), m_Ny(ny), m_grid(nx * ny) {}

  size_t index(size_t ix, size_t iy) const { return iy + m_Ny * ix; }

  T& get(size_t ix, size_t iy) { return m_grid[index(ix, iy)]; }
  const T& get(size_t ix, size_t iy) const { return m_grid[index(ix, iy)]; }

  std::vector<T>& get() { return m_grid; }
  const std::vector<T>& get() const { return m_grid; }

  T max() const { return *std::max_element(m_grid.begin(), m_grid.end()); }
  T min() const { return *std::min_element(m_grid.begin(), m_grid.end()); }

  size_t getNx() const { return m_Nx; }
  size_t getNy() const { return m_Ny; }
  size_t size() const { return m_Nx * m_Ny; }

  void copy(std::vector<T> v) {
    assert(v.size() == m_grid.size());
    m_grid = std::move(v);
  }

  template <typename Func>
  void for_each(Func fn) {
    std::for_each(m_grid.begin(), m_grid.end(), fn);
  }

 private:
  size_t m_Nx = 0;
  size_t m_Ny = 0;
  std::vector<T> m_grid;
};

}  // namespace CRAMS

#endif  // CRAMS_UTILS_GRID_H_
