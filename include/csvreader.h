#ifndef INCLUDE_CSVREADER_H__
#define INCLUDE_CSVREADER_H__

#include <fstream>
#include <string>
#include <vector>

namespace CRAMS {

class CSVReader {
  std::string fileName;
  std::string delimeter;

 public:
  CSVReader(std::string filename, std::string delm = ",") : fileName(filename), delimeter(delm) {}
  std::vector<std::vector<std::string> > getData() const;
};

}  // namespace CRAMS

#endif