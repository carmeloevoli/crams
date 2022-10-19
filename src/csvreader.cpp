#include "csvreader.h"

#include "utilities.h"

namespace CRAMS {

std::vector<std::string> split(const std::string& str, const std::string& delim) {
  if (str == "" || delim == "") throw std::invalid_argument("str or delim cannot be empty");

  std::vector<std::string> tokens;
  size_t prev = 0, pos = 0;
  do {
    pos = str.find(delim, prev);
    if (pos == std::string::npos) pos = str.length();
    auto token = str.substr(prev, pos - prev);
    if (!token.empty()) tokens.push_back(token);
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());
  return tokens;
}

std::vector<std::vector<std::string> > CSVReader::getData() const {
  std::ifstream file(fileName);
  std::vector<std::vector<std::string> > dataList;
  std::string line = "";
  const char vetoChar = '#';
  while (getline(file, line)) {
    if (line.at(0) != vetoChar) {
      auto vec = split(line, delimeter);
      dataList.push_back(vec);
    }
  }
  file.close();
  return dataList;
}

}  // namespace CRAMS