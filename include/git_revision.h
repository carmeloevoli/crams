#ifndef INCLUDE_GIT_REVISION_H_
#define INCLUDE_GIT_REVISION_H_

#include <string>

std::string git_sha1();
std::string get_version();
bool git_has_local_changes();


#endif /* INCLUDE_GIT_REVISION_H_ */
