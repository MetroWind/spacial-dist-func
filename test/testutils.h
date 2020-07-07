// -*- mode: c++; -*-
#ifndef SDF_TEST_UTILS_H
#define SDF_TEST_UTILS_H

#include <random>

#include <Eigen/Dense>

namespace TestGlobal
{
    static std::random_device Dev;
    static std::mt19937 RandAlgo(Dev());
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
