#ifndef CRAMS_UTILS_GIT_REVISION_H_
#define CRAMS_UTILS_GIT_REVISION_H_

#include <string>

namespace CRAMS {

std::string git_sha1();
std::string get_version();
bool git_has_local_changes();

}  // namespace CRAMS

#endif  // CRAMS_UTILS_GIT_REVISION_H_
