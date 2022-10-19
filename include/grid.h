#ifndef INCLUDE_XS4GCR_GRID_H
#define INCLUDE_XS4GCR_GRID_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <vector>

namespace CRAMS {

/** Lower and upper neighbor in a unit grid */
inline void clamp(double x, int n, int& lo, int& hi) {
  lo = int(floor(x * (n - 1)));
  hi = (lo + 1) % n;
}

template <typename T>
class Grid {
 protected:
  std::vector<T> m_grid;
  size_t m_Nx = 0;
  size_t m_Ny = 0;

 public:
  Grid() {}

  Grid(size_t nx, size_t ny) : m_Nx(nx), m_Ny(ny) { m_grid.resize(nx * ny); }

  virtual ~Grid() = default;

  size_t index(size_t ix, size_t iy) const { return iy + m_Ny * ix; }

  /** Accessor / Mutator */
  T& get(size_t ix, size_t iy) { return m_grid[index(ix, iy)]; }

  T& get(size_t i) { return m_grid[i]; }

  /* Min / Max */
  T max() const { return *max_element(m_grid.begin(), m_grid.end()); }

  T min() const { return *min_element(m_grid.begin(), m_grid.end()); }

  /** Accessor */
  const T& get(size_t ix, size_t iy) const { return m_grid[index(ix, iy)]; }

  const T& get(size_t i) const { return m_grid[i]; }

  /** Return a reference to the grid values */
  std::vector<T>& get() { return m_grid; }

  const std::vector<T>& get() const { return m_grid; }

  void clear() { m_grid.clear(); }

  size_t getNX() const { return m_Nx; }

  size_t getNY() const { return m_Ny; }

  size_t size() const { return (m_Nx * m_Ny); }

  /** Interpolate the grid at a given position */
  T interpolate(double x, double y) const {
    // position on a unit grid
    // Vector3d r = (position - gridOrigin) / spacing;

    // indices of lower and upper neighbors
    int ix, iX, iy, iY;
    clamp(x, m_Nx, ix, iX);
    clamp(y, m_Ny, iy, iY);

    // linear fraction to lower and upper neighbors
    double fx = x - floor(x);
    double fX = 1 - fx;
    double fy = y - floor(y);
    double fY = 1 - fy;

    // bilinear interpolation (see http://paulbourke.net/miscellaneous/interpolation)
    T b(0.);
    // V00 (1 - x) (1 - y) +
    b += get(ix, iy) * fX * fY;
    // V10 x (1 - y) +
    b += get(iX, iy) * fx * fY;
    // V01 (1 - x) y +
    b += get(ix, iY) * fX * fy;
    // V11 x y
    b += get(iX, iY) * fx * fy;

    return b;
  }

  void copy(std::vector<T> v) {
    assert(v.size() == m_grid.size());
    m_grid.assign(v.begin(), v.end());
  }

  void for_each(std::function<void(T&)> fn) { std::for_each(m_grid.begin(), m_grid.end(), fn); }
};

}  // namespace CRAMS

#endif  // INCLUDE_XS4GCR_GRID_H