#include "fs.hxx"
#include "../base.hxx"

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/unistd.h>
#include <dirent.h>
#include <string.h>
#include <sstream>

namespace desres { namespace viparr { namespace fs {

    bool exists(std::string const& path) {
        struct stat buf[1];
        return stat(path.data(), buf) == 0;
    }

    bool is_directory(std::string const& path) {
        struct stat buf[1];
        return stat(path.data(), buf)==0 && S_ISDIR(buf->st_mode);
    }

    bool is_empty(std::string const& path) {
        struct dirent *d;
        DIR *dir = opendir(path.data());
        if (!dir) return true;
        int n = 0;
        while ((d = readdir(dir)) != NULL) {
            if (++n > 2) break;
        }
        closedir(dir);
        return n<=2;
    }

    std::vector<std::string> iter_directory(std::string const& path) {
        struct dirent *d;
        std::vector<std::string> entries;
        DIR *dir = opendir(path.data());
        if (dir) {
            while ((d = readdir(dir)) != NULL) {
                const char* name = d->d_name;
                if (name[0]=='.' && (name[1] == '\0' || name[1] == '.')) continue;
                entries.push_back(name);
            }
        }
        closedir(dir);
        return entries;
    }

    void create_directories(std::string const& dir) {
        mode_t mode = 0755;
        /* Adapted from http://stackoverflow.com/a/2336245/119527 */
        const char* path = dir.data();
        const size_t len = dir.size();
        char _path[PATH_MAX];
        char *p;
    
        errno = 0;
    
        /* Copy string so its mutable */
        if (len > sizeof(_path)-1) {
            errno = ENAMETOOLONG;
            VIPARR_FAIL(dir << " : " << strerror(errno));
        }
        strcpy(_path, path);
    
        /* Iterate the string */
        for (p = _path + 1; *p; p++) {
            if (*p == '/') {
                /* Temporarily truncate */
                *p = '\0';
    
                if (mkdir(_path, mode) != 0) {
                    if (errno != EEXIST)
                        VIPARR_FAIL(dir << " : " << strerror(errno));
                }
    
                *p = '/';
            }
        }
    
        if (mkdir(_path, mode) != 0) {
            if (errno != EEXIST)
                VIPARR_FAIL(dir << " : " << strerror(errno));
        }
    }


}}}

