#include "crams/utils/csvreader.h"

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::vector<std::string> split(const std::string& str, const std::string& delim) {
  if (str.empty() || delim.empty())
    throw std::invalid_argument("CSVReader: str and delim must be non-empty");
  std::vector<std::string> tokens;
  size_t prev = 0, pos = 0;
  do {
    pos = str.find(delim, prev);
    if (pos == std::string::npos) pos = str.length();
    auto token = str.substr(prev, pos - prev);
    if (!token.empty()) tokens.push_back(std::move(token));
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());
  return tokens;
}

}  // namespace

namespace CRAMS {

CSVReader::CSVReader(std::string filename, std::string delimiter)
    : m_filename(std::move(filename)), m_delimiter(std::move(delimiter)) {}

std::vector<std::vector<std::string>> CSVReader::getData() const {
  std::ifstream file(m_filename);
  if (!file.is_open())
    throw std::runtime_error("CSVReader: cannot open file '" + m_filename + "'");

  std::vector<std::vector<std::string>> dataList;
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    dataList.push_back(split(line, m_delimiter));
  }
  return dataList;
}

}  // namespace CRAMS
