#include "crams/fragmentation.h"

#include <plog/Log.h>

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "crams/utils/csvreader.h"
#include "crams/utils/numeric.h"
#include "crams/utils/utilities.h"

namespace CRAMS {

namespace {

// Checks that the table's first-row energy grid (in GeV, after *idColumns* label
// cells) matches the grid the reader built from T_min/T_max/T_size.
void validateEnergyGrid(const std::string& modelName, const std::vector<std::string>& header, size_t idColumns,
                        const std::vector<double>& grid) {
  if (header.size() != idColumns + grid.size())
    throw std::runtime_error(modelName + ": energy grid header has " + std::to_string(header.size()) +
                             " columns, expected " + std::to_string(idColumns + grid.size()));
  constexpr double tol = 1e-4;
  for (size_t i = 0; i < grid.size(); ++i) {
    double fileT;
    try {
      fileT = std::stod(header[idColumns + i]);
    } catch (const std::exception&) {
      throw std::runtime_error(modelName + ": non-numeric energy grid value '" + header[idColumns + i] + "'");
    }
    const double codeT = grid[i] / CGS::GeV;
    if (std::abs(fileT - codeT) > tol * codeT)
      throw std::runtime_error(modelName + ": energy grid mismatch at column " + std::to_string(idColumns + i) +
                               " (file " + std::to_string(fileT) + " GeV vs built " + std::to_string(codeT) + " GeV)");
  }
}

}  // namespace

NucFragXsec::~NucFragXsec() { LOGD << "deleted NucFragXsec"; }

double NucFragXsec::getXsecOnISM(const PID& projectile, const PID& fragment, const double& T) const {
  const double sigma_H = getXsecOnHtarget(projectile, fragment, T);
  return sigma_H * (1. + CGS::K_He * CGS::f_He) / (1. + CGS::f_He);
}

NucFragFromTable::NucFragFromTable(std::string modelName, std::string tableFilename, double T_min, double T_max,
                                   size_t T_size)
    : m_modelName(std::move(modelName)),
      m_tableFilename(std::move(tableFilename)),
      m_T_min(T_min),
      m_T_max(T_max),
      m_T_size(T_size) {
  buildEnergyArray();
  if (Utilities::fileExists(m_tableFilename))
    loadXsecTable();
  else
    throw std::runtime_error(m_modelName + " fragmentation xsecs file not found: " + m_tableFilename);
  if (m_table.empty())
    throw std::runtime_error(m_modelName + " fragmentation xsecs table is empty: " + m_tableFilename);
  LOGD << m_modelName << " fragmentation table read with " << m_table.size() << " channels";
}

void NucFragFromTable::buildEnergyArray() {
  const double logRatio = std::log(m_T_max / m_T_min) / (m_T_size - 1);
  for (size_t i = 0; i < m_T_size; ++i) m_T.push_back(m_T_min * std::exp(i * logRatio));
}

double NucFragFromTable::getXsecOnHtarget(const PID& projectile, const PID& fragment, const double& T) const {
  const auto it = m_table.find({projectile, fragment});
  if (it == m_table.end()) return 0.;  // channel not in this model -> no contribution
  // Small relative tolerance at the edges: the simulation energy grid can land
  // a hair outside [m_T_min, m_T_max] due to floating-point rounding in the
  // log-spaced grids, which is not a genuine out-of-range request.
  constexpr double edgeTol = 1e-6;
  if (T < m_T_min * (1. - edgeTol) || T > m_T_max * (1. + edgeTol))
    throw std::runtime_error(m_modelName + " fragmentation xsec requested at T = " + std::to_string(T / CGS::GeV) +
                             " GeV, outside the tabulated range [" + std::to_string(m_T_min / CGS::GeV) + ", " +
                             std::to_string(m_T_max / CGS::GeV) + "] GeV");
  // Clamp into the tabulated grid so a boundary point does not trip the interpolator.
  const double T_clamped = std::min(std::max(T, m_T.front()), m_T.back());
  return Numeric::LinearInterpolator<double>(m_T, it->second, T_clamped);
}

void NucFragFromTable::loadXsecTable() {
  const auto headerAndData = CSVReader(m_tableFilename).getHeaderAndData();
  validateEnergyGrid(m_modelName, headerAndData.first, 4, m_T);  // Z_proj, A_proj, Z_frag, A_frag + energy grid

  const size_t expectedColumns = m_T_size + 4;  // Z_proj, A_proj, Z_frag, A_frag, sigma(T_0 .. T_N-1)
  for (const auto& row : headerAndData.second) {
    if (row.size() != expectedColumns)
      throw std::runtime_error("malformed " + m_modelName + " fragmentation xsecs row: expected " +
                               std::to_string(expectedColumns) + " columns, got " + std::to_string(row.size()));
    const PID projectile(static_cast<int>(row[0]), static_cast<int>(row[1]));
    const PID fragment(static_cast<int>(row[2]), static_cast<int>(row[3]));
    std::vector<double> x;
    x.reserve(m_T_size);
    for (size_t i = 0; i < m_T_size; ++i) x.push_back(row[4 + i] * CGS::mbarn);
    m_table[{projectile, fragment}] = x;
  }
}

constexpr double kFragTmin = 0.01 * CGS::GeV;
constexpr double kFragTmax = 1e5 * CGS::GeV;
constexpr size_t kFragTsize = 112;

NucFragFluka4Dragon::NucFragFluka4Dragon()
    : NucFragFromTable("Fluka4Dragon", "data/crams_fragmentation_fluka4dragon.csv", kFragTmin, kFragTmax, kFragTsize) {}

NucFragUsineGalprop17Opt12::NucFragUsineGalprop17Opt12()
    : NucFragFromTable("USINE_GALPROP17_OPT12", "data/crams_fragmentation_usine_galprop17_opt12.csv", kFragTmin,
                       kFragTmax, kFragTsize) {}

NucFragUsineGalprop17Opt22::NucFragUsineGalprop17Opt22()
    : NucFragFromTable("USINE_GALPROP17_OPT22", "data/crams_fragmentation_usine_galprop17_opt22.csv", kFragTmin,
                       kFragTmax, kFragTsize) {}

NucFragUsineWebber03Coste12::NucFragUsineWebber03Coste12()
    : NucFragFromTable("USINE_WEBBER03_COSTE12", "data/crams_fragmentation_usine_webber03+coste12.csv", kFragTmin,
                       kFragTmax, kFragTsize) {}

}  // namespace CRAMS
