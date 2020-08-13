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

#include "sdf.h"
#include "utils.h"

using namespace libmd;

namespace libmd
{
    Eigen::Matrix3f rotateToAlignX(const VecRefType& vec,
                                   const VecRefType& to_xy)
    {
        Eigen::Matrix3f Rot1;
        {
            Eigen::Vector3f Axis;
            Axis[0] = 0.0f;
            Axis[1] = vec[2];
            Axis[2] = -vec[1];
            Axis /= Axis.norm();

            const float Cos = vec[0] / vec.norm();
            const float Sin = std::sqrt(1.0f - Cos*Cos);

            Rot1(0, 0) = Cos;
            Rot1(0, 1) = -Axis[2] * Sin;
            Rot1(0, 2) = Axis[1] * Sin;
            Rot1(1, 0) = Axis[2] * Sin;
            Rot1(1, 1) = Cos + Axis[1] * Axis[1] * (1.0f - Cos);
            Rot1(1, 2) = Axis[1] * Axis[2] * (1.0f - Cos);
            Rot1(2, 0) = -Axis[1] * Sin;
            Rot1(2, 1) = Axis[2] * Axis[1] * (1.0f - Cos);
            Rot1(2, 2) = Cos + Axis[2] * Axis[2] * (1.0f - Cos);
        }

        const Eigen::Vector3f Rotated = Rot1 * to_xy;
        Eigen::Matrix3f Rot2 = Eigen::Matrix3f::Zero();
        {
            float Cos = Rotated[1] / std::sqrt(Rotated[1] * Rotated[1] +
                                                     Rotated[2] * Rotated[2]);
            float Sin = std::sqrt(1.0f - Cos*Cos);
            if(Rotated[2] > 0)
            {
                Sin = -Sin;
            }

            Rot2(0, 0) = 1.0f;
            Rot2(1, 1) = Cos;
            Rot2(1, 2) = -Sin;
            Rot2(2, 1) = Sin;
            Rot2(2, 2) = Cos;
        }

        return Rot2 * Rot1;
    }
} // namespace libmd

namespace sdf
{
    const typename DistCountTraits::ValueType DistCountTraits::Zero = 0;
    const typename DistChargeTraits::ValueType DistChargeTraits::Zero = 0;
    const typename DistDetailedCountTraits::ValueType DistDetailedCountTraits::Zero = {};
} // namespace sdf
