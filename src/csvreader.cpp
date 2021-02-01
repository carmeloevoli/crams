#include "csvreader.h"

#include "utilities.h"

namespace CRAMS {

std::vector<std::vector<std::string> > CSVReader::getData() const {
  std::ifstream file(fileName);
  std::vector<std::vector<std::string> > dataList;
  std::string line = "";
  const char vetoChar = '#';
  while (getline(file, line)) {
    if (line.at(0) != vetoChar) {
      auto vec = Utilities::split(line, delimeter);
      dataList.push_back(vec);
    }
  }
  file.close();
  return dataList;
}

}  // namespace CRAMS