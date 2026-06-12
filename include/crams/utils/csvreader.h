#ifndef CRAMS_UTILS_CSVREADER_H_
#define CRAMS_UTILS_CSVREADER_H_

#include <string>
#include <utility>
#include <vector>

namespace CRAMS {

class CSVReader {
 public:
  explicit CSVReader(std::string filename, std::string delimiter = ",");
  std::vector<std::vector<std::string>> getData() const;

  // Like getData() but converts every cell to double as it streams, without
  // holding the whole file as strings. Throws on any non-numeric cell.
  std::vector<std::vector<double>> getDataAsDouble() const;

  // For tables whose first non-comment row holds column labels/coordinates and
  // the rest are numeric: returns {header (as strings), data (as doubles)} in a
  // single streaming pass. Throws on a non-numeric data cell or if no rows.
  std::pair<std::vector<std::string>, std::vector<std::vector<double>>> getHeaderAndData() const;

 private:
  std::string m_filename;
  std::string m_delimiter;
};

}  // namespace CRAMS

#endif  // CRAMS_UTILS_CSVREADER_H_
