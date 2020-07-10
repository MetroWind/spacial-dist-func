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

#include <cmath>

#include "pbc.h"

namespace libmd
{
    float RectPbc3d :: dist1d(const size_t dim_idx, float lhs, float rhs) const
    {
        float Dim = Dimension[dim_idx];
        float Dist = std::abs(lhs - rhs);
        if(Dist > Dim)
        {
            Dist -= std::floor(Dist / Dim) * Dim;
        }

        if(Dist > Dim * 0.5)
        {
            return Dim - Dist;
        }
        else
        {
            return Dist;
        }
    }

    float RectPbc3d :: distSquare(const VecRefType& lhs,
                                  const VecRefType& rhs) const
    {
        float x = dist1d(0, lhs[0], rhs[0]);
        float y = dist1d(1, lhs[1], rhs[1]);
        float z = dist1d(2, lhs[2], rhs[2]);

        return x*x + y*y + z*z;
    }

    float RectPbc3d :: wrap1d(const size_t dim_idx, float base, float rhs) const
    {
        const float Dim = Dimension[dim_idx];
        float Dist = rhs - base;
        if(Dist > Dim || -Dist >= Dim)
        {
            rhs -= std::floor(Dist / Dim) * Dim;
        }
        if(rhs - base > Dim * 0.5)
        {
            rhs -= Dim;
        }
        else if(base - rhs > Dim * 0.5) // See test case “PBC wrap forward”.
        {
            rhs += Dim;
        }

        return rhs;
    }

    void RectPbc3d :: wrapVec(const float base[], float to_wrap[]) const
    {
        to_wrap[0] = wrap1d(0, base[0], to_wrap[0]);
        to_wrap[1] = wrap1d(1, base[1], to_wrap[1]);
        to_wrap[2] = wrap1d(2, base[2], to_wrap[2]);
    }

} // namespace libmd
