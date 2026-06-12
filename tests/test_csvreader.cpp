#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "crams/utils/csvreader.h"

static int g_pass = 0;
static int g_fail = 0;

#define CHECK(cond)                                                                    \
  do {                                                                                 \
    if (cond) {                                                                        \
      ++g_pass;                                                                        \
    } else {                                                                           \
      ++g_fail;                                                                        \
      std::cerr << "FAIL: " << #cond << " at " << __FILE__ << ":" << __LINE__ << "\n"; \
    }                                                                                  \
  } while (0)

#define CHECK_THROW(expr, exc) \
  do {                         \
    bool caught = false;       \
    try {                      \
      (void)(expr);            \
    } catch (const exc&) {     \
      caught = true;           \
    }                          \
    CHECK(caught);             \
  } while (0)

static void write_file(const std::string& path, const std::string& content) {
  std::ofstream f(path);
  assert(f.is_open());
  f << content;
}

static const std::string TMP = "/tmp/test_csvreader_";

void test_basic_csv() {
  const std::string path = TMP + "basic.csv";
  write_file(path, "1,2,3\n4,5,6\n");
  CRAMS::CSVReader reader(path);
  auto data = reader.getData();
  CHECK(data.size() == 2);
  CHECK(data[0].size() == 3);
  CHECK(data[0][0] == "1");
  CHECK(data[0][1] == "2");
  CHECK(data[0][2] == "3");
  CHECK(data[1][0] == "4");
  CHECK(data[1][1] == "5");
  CHECK(data[1][2] == "6");
}

void test_comment_lines_skipped() {
  const std::string path = TMP + "comments.csv";
  write_file(path, "# this is a comment\n1,2,3\n# another comment\n4,5,6\n");
  CRAMS::CSVReader reader(path);
  auto data = reader.getData();
  CHECK(data.size() == 2);
  CHECK(data[0][0] == "1");
  CHECK(data[1][0] == "4");
}

void test_empty_lines_skipped() {
  const std::string path = TMP + "emptylines.csv";
  write_file(path, "\n1,2,3\n\n4,5,6\n\n");
  CRAMS::CSVReader reader(path);
  auto data = reader.getData();
  CHECK(data.size() == 2);
  CHECK(data[0][0] == "1");
  CHECK(data[1][0] == "4");
}

void test_mixed_comments_and_empty_lines() {
  const std::string path = TMP + "mixed.csv";
  write_file(path, "# header\n\n1,2\n# comment\n\n3,4\n");
  CRAMS::CSVReader reader(path);
  auto data = reader.getData();
  CHECK(data.size() == 2);
  CHECK(data[0][0] == "1");
  CHECK(data[1][0] == "3");
}

void test_custom_delimiter_space() {
  const std::string path = TMP + "space.csv";
  write_file(path, "1 2 3\n4 5 6\n");
  CRAMS::CSVReader reader(path, " ");
  auto data = reader.getData();
  CHECK(data.size() == 2);
  CHECK(data[0][0] == "1");
  CHECK(data[0][1] == "2");
  CHECK(data[0][2] == "3");
}

void test_custom_delimiter_tab() {
  const std::string path = TMP + "tab.csv";
  write_file(path, "a\tb\tc\n");
  CRAMS::CSVReader reader(path, "\t");
  auto data = reader.getData();
  CHECK(data.size() == 1);
  CHECK(data[0].size() == 3);
  CHECK(data[0][0] == "a");
  CHECK(data[0][1] == "b");
  CHECK(data[0][2] == "c");
}

void test_single_column() {
  const std::string path = TMP + "single.csv";
  write_file(path, "hello\nworld\n");
  CRAMS::CSVReader reader(path);
  auto data = reader.getData();
  CHECK(data.size() == 2);
  CHECK(data[0][0] == "hello");
  CHECK(data[1][0] == "world");
}

void test_file_not_found_throws() {
  CRAMS::CSVReader reader("/tmp/this_file_does_not_exist_xyzzy.csv");
  CHECK_THROW(reader.getData(), std::runtime_error);
}

void test_only_comments_returns_empty() {
  const std::string path = TMP + "onlycomments.csv";
  write_file(path, "# line 1\n# line 2\n");
  CRAMS::CSVReader reader(path);
  auto data = reader.getData();
  CHECK(data.empty());
}

void test_no_trailing_newline() {
  const std::string path = TMP + "nonewline.csv";
  write_file(path, "1,2,3");
  CRAMS::CSVReader reader(path);
  auto data = reader.getData();
  CHECK(data.size() == 1);
  CHECK(data[0].size() == 3);
}

void test_numeric_values() {
  const std::string path = TMP + "numeric.csv";
  write_file(path, "1,3.14,-2.71\n0,1e10,2.5e-3\n");
  CRAMS::CSVReader reader(path);
  auto data = reader.getData();
  CHECK(data.size() == 2);
  CHECK(std::stod(data[0][1]) == 3.14);
  CHECK(std::stod(data[1][2]) == 2.5e-3);
}

void test_getDataAsDouble() {
  const std::string path = TMP + "asdouble.csv";
  write_file(path, "# header\n1,3.14,-2.71\n0,1e10,2.5e-3\n");
  CRAMS::CSVReader reader(path);
  auto data = reader.getDataAsDouble();
  CHECK(data.size() == 2);
  CHECK(data[0].size() == 3);
  CHECK(data[0][0] == 1.0);
  CHECK(data[0][1] == 3.14);
  CHECK(data[1][1] == 1e10);
  CHECK(data[1][2] == 2.5e-3);
}

void test_getDataAsDouble_tab_delimited() {
  const std::string path = TMP + "asdouble_tab.csv";
  write_file(path, "1\t2\t3.5\n");
  CRAMS::CSVReader reader(path, "\t");
  auto data = reader.getDataAsDouble();
  CHECK(data.size() == 1);
  CHECK(data[0].size() == 3);
  CHECK(data[0][2] == 3.5);
}

void test_getDataAsDouble_non_numeric_throws() {
  const std::string path = TMP + "asdouble_bad.csv";
  write_file(path, "1,2,3\n4,oops,6\n");
  CRAMS::CSVReader reader(path);
  CHECK_THROW(reader.getDataAsDouble(), std::runtime_error);
}

void test_getHeaderAndData() {
  const std::string path = TMP + "headerdata.csv";
  write_file(path, "# comment\nZ,A,0.01,0.02\n1,1,3.5,4.5\n2,4,6.0,7.0\n");
  CRAMS::CSVReader reader(path);
  auto result = reader.getHeaderAndData();
  const auto& header = result.first;
  const auto& data = result.second;
  CHECK(header.size() == 4);
  CHECK(header[0] == "Z");
  CHECK(std::stod(header[2]) == 0.01);
  CHECK(data.size() == 2);
  CHECK(data[0].size() == 4);
  CHECK(data[1][2] == 6.0);
}

void test_getHeaderAndData_no_rows_throws() {
  const std::string path = TMP + "headeronly_empty.csv";
  write_file(path, "# only comments\n");
  CRAMS::CSVReader reader(path);
  CHECK_THROW(reader.getHeaderAndData(), std::runtime_error);
}

int main() {
  test_basic_csv();
  test_comment_lines_skipped();
  test_empty_lines_skipped();
  test_mixed_comments_and_empty_lines();
  test_custom_delimiter_space();
  test_custom_delimiter_tab();
  test_single_column();
  test_file_not_found_throws();
  test_only_comments_returns_empty();
  test_no_trailing_newline();
  test_numeric_values();
  test_getDataAsDouble();
  test_getDataAsDouble_tab_delimited();
  test_getDataAsDouble_non_numeric_throws();
  test_getHeaderAndData();
  test_getHeaderAndData_no_rows_throws();

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
