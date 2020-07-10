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

#ifndef SDF_ENDIAN_H
#define SDF_ENDIAN_H

#include <algorithm>

namespace libmd
{
    namespace Endian
    {
        enum Endian { BIG, LITTLE };
        inline Endian current()
        {
            union{ uint16_t i; char c[2]; } bint = {0x0102};

            if(bint.c[0] == 1)
            {
                return BIG;
            }
            else
            {
                return LITTLE;
            }
        }

        template <typename T> T swapped(const T& x)
        {
            union U
            {
                T val;
                std::array<std::uint8_t, sizeof(T)> raw;
            } src, dst;

            src.val = x;
            std::reverse_copy(src.raw.begin(), src.raw.end(), dst.raw.begin());
            return dst.val;
        }

        template <typename T> void swap(T& x)
        {
            union U
            {
                T val;
                std::array<std::uint8_t, sizeof(T)> raw;
            } src, dst;
            src.val = x;
            std::reverse_copy(src.raw.begin(), src.raw.end(), dst.raw.begin());
            x = dst.val;
        }
    }
}

#endif
