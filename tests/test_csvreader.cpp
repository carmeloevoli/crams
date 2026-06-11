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

  std::cout << g_pass << " passed, " << g_fail << " failed\n";
  return g_fail > 0 ? 1 : 0;
}
