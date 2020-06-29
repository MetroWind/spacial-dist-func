#ifndef SDF_PBC_H
#define SDF_PBC_H

#include <iostream>
#include <array>
#include <cmath>

#include <Eigen/Dense>

namespace libmd
{
    class RectPbc3d
    {
    public:
        RectPbc3d() = delete;
        RectPbc3d(float x, float y, float z) : Dimension({ x, y, z }) {}

        using VecRefType = Eigen::Ref<Eigen::Vector3f>;
        float dist(const VecRefType& lhs, const VecRefType& rhs) const
        {
            return std::sqrt(distSquare(lhs, rhs));
        }
        float distSquare(const VecRefType& lhs, const VecRefType& rhs) const;

    private:
        float dist1d(const size_t dim, float lhs, float rhs) const;
        const std::array<float, 3> Dimension;

    };

} // namespace libmd

#endif
