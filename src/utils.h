// -*- mode: c++; -*-
#ifndef SDF_UTILS_H
#define SDF_UTILS_H

#include <Eigen/Dense>

#define UNUSED(x) (void)(x)

namespace libmd
{
    using V3Map = Eigen::Map<Eigen::Vector3f>;
    using VecRefType = Eigen::Ref<Eigen::Vector3f>;
}

#endif
