#ifndef viparr_importexport_fs
#define viparr_importexport_fs

#include <string>
#include <vector>

namespace desres { namespace viparr { namespace fs {

    // filesystem wrappers
    bool exists(std::string const& path);
    bool is_directory(std::string const& path);
    bool is_empty(std::string const& dir);
    std::vector<std::string> iter_directory(std::string const& dir);
    void create_directories(std::string const& dir);    // throws

}}}

#endif


