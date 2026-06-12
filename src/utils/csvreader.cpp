#include "crams/utils/csvreader.h"

#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

std::vector<std::string> split(const std::string& str, const std::string& delim) {
  if (str.empty() || delim.empty()) throw std::invalid_argument("CSVReader: str and delim must be non-empty");
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
  if (!file.is_open()) throw std::runtime_error("CSVReader: cannot open file '" + m_filename + "'");

  std::vector<std::vector<std::string>> dataList;
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    dataList.push_back(split(line, m_delimiter));
  }
  return dataList;
}

namespace {

std::vector<double> toDoubleRow(const std::vector<std::string>& tokens, const std::string& filename, size_t rowIndex) {
  std::vector<double> row;
  row.reserve(tokens.size());
  for (const auto& token : tokens) {
    try {
      row.push_back(std::stod(token));
    } catch (const std::exception&) {
      throw std::runtime_error("CSVReader: non-numeric value '" + token + "' in '" + filename + "' data row " +
                               std::to_string(rowIndex));
    }
  }
  return row;
}

}  // namespace

std::vector<std::vector<double>> CSVReader::getDataAsDouble() const {
  std::ifstream file(m_filename);
  if (!file.is_open()) throw std::runtime_error("CSVReader: cannot open file '" + m_filename + "'");

  std::vector<std::vector<double>> dataList;
  std::string line;
  size_t rowIndex = 0;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    dataList.push_back(toDoubleRow(split(line, m_delimiter), m_filename, ++rowIndex));
  }
  return dataList;
}

std::pair<std::vector<std::string>, std::vector<std::vector<double>>> CSVReader::getHeaderAndData() const {
  std::ifstream file(m_filename);
  if (!file.is_open()) throw std::runtime_error("CSVReader: cannot open file '" + m_filename + "'");

  std::vector<std::string> header;
  std::vector<std::vector<double>> dataList;
  std::string line;
  bool haveHeader = false;
  size_t rowIndex = 0;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (!haveHeader) {
      header = split(line, m_delimiter);
      haveHeader = true;
      continue;
    }
    dataList.push_back(toDoubleRow(split(line, m_delimiter), m_filename, ++rowIndex));
  }
  if (!haveHeader) throw std::runtime_error("CSVReader: no rows in '" + m_filename + "'");
  return {std::move(header), std::move(dataList)};
}

}  // namespace CRAMS
