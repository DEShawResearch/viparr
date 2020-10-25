#ifndef viparr_util_util_h
#define viparr_util_util_h

#include <vector>
#include <string>

namespace desres { namespace viparr {

    std::vector<std::string> ViparrSplitString(std::string word, char sep = ' ');
    void ViparrReplaceAll(std::string& word, std::string const& oldstr, std::string const& newstr);

}}

#endif
