#include "utils.h"

namespace libmd
{
    // Equivalent of Python’s str.strip().
    std::string strip(const std::string& str,
                      const std::string& whitespace /* = " \t\n" */)
    {
        const auto StrBegin = str.find_first_not_of(whitespace);
        if (StrBegin == std::string::npos)
            return ""; // no content

        const auto StrEnd = str.find_last_not_of(whitespace);
        const auto StrRange = StrEnd - StrBegin + 1;

        return str.substr(StrBegin, StrRange);
    }
}
