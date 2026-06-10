#ifndef CRAMS_UTILS_CSVREADER_H_
#define CRAMS_UTILS_CSVREADER_H_

#include <string>
#include <vector>

namespace CRAMS {

class CSVReader {
 public:
  explicit CSVReader(std::string filename, std::string delimiter = ",");
  std::vector<std::vector<std::string>> getData() const;

 private:
  std::string m_filename;
  std::string m_delimiter;
};

}  // namespace CRAMS

#endif  // CRAMS_UTILS_CSVREADER_H_
