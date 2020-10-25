#include "util.hxx"
#include <sstream>
#include <iterator>
#include <algorithm>

std::vector<std::string> 
desres::viparr::ViparrSplitString(std::string sentence, char sep) {
    std::istringstream iss(sentence);
    std::vector<std::string> words;
    std::copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter(words));
    return words;
}

void
desres::viparr::ViparrReplaceAll(std::string& word, std::string const& oldstr, std::string const& newstr) {
    size_t index = 0;
    for (;;) {
        index = word.find(oldstr, index);
        if (index == std::string::npos) break;
        word.replace(index, oldstr.size(), newstr);
        index += newstr.size();
    }
}

