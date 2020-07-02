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
        using VecRefType = Eigen::Ref<Eigen::Vector3f>;

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
