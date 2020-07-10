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

#ifndef SDF_PBC_H
#define SDF_PBC_H

#include <iostream>
#include <array>
#include <cmath>

#include <Eigen/Dense>

#include "utils.h"

namespace libmd
{
    class RectPbc3d
    {
    public:
        RectPbc3d() = delete;
        RectPbc3d(float x, float y, float z)
                : DiagLength(std::sqrt(x*x + y*y + z*z)),
                  Dimension({ x, y, z }) {}

        float dist(const VecRefType& lhs, const VecRefType& rhs) const
        {
            return std::sqrt(distSquare(lhs, rhs));
        }
        float distSquare(const VecRefType& lhs, const VecRefType& rhs) const;

        void wrapVec(const float base[], float to_wrap[]) const;

        // This really should be
        //
        //   void wrapVec(const VecRefType&, VecRefType) const
        //
        // But for some reason it does not compile with both arguments
        // being Vector3f: “calling a private constructor of class
        // 'Eigen::Ref<Eigen::Matrix<float, 3, 1, 0, 3, 1>, 0,
        // Eigen::InnerStride<1> >'”...
        template <typename T>
        void wrapVec(const T& base, T& to_wrap) const
        {
            wrapVec(base.data(), to_wrap.data());
        }

        const float DiagLength;

    private:
        float dist1d(const size_t dim_idx, float lhs, float rhs) const;
        float wrap1d(const size_t dim_idx, float base, float rhs) const;

        const std::array<float, 3> Dimension;
    };

} // namespace libmd

#endif
