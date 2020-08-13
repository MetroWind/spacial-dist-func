// -*- mode: c++; -*-
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

#ifndef SDF_TEST_UTILS_H
#define SDF_TEST_UTILS_H

#include <random>

#include <Eigen/Dense>

namespace TestGlobal
{
    static std::random_device Dev;
    static std::mt19937 RandAlgo(Dev());
    constexpr float FloatMargin = 0.0002;
}

inline float randUni(float low, float high)
{
    std::uniform_real_distribution<float> Dist(low, high);
    return Dist(TestGlobal::RandAlgo);
}

inline Eigen::Vector3f randVec(std::pair<float, float> rangex,
                        std::pair<float, float> rangey,
                        std::pair<float, float> rangez)
{
    return Eigen::Vector3f(randUni(rangex.first, rangex.second),
                           randUni(rangey.first, rangey.second),
                           randUni(rangez.first, rangez.second));
}

#endif
