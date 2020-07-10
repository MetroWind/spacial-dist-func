// Copyright 2020 MetroWind <chris.corsair@gmail.com>
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see
// <https://www.gnu.org/licenses/>.

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
