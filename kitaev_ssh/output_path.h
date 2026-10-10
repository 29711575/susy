#ifndef KITAEV_SSH_OUTPUT_PATH_H
#define KITAEV_SSH_OUTPUT_PATH_H
#include <cerrno>
#include <string>
#include <system_error>
#include <sys/stat.h>
namespace kitaev_ssh {
// Linux HPC builds may have a C++17 compiler but pre-filesystem libstdc++.
// lstat also protects dangling symlinks from accidental output replacement.
inline bool outputPathExists(const std::string &path) {
    struct stat entry;
    if (::lstat(path.c_str(), &entry)==0) return true;
    const int error=errno;
    if (error==ENOENT) return false;
    throw std::system_error(error,std::generic_category(),"inspect output path: "+path);
}
}
#endif
