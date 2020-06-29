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

    float RectPbc3d :: distSquare(const RectPbc3d::VecRefType& lhs,
                                  const RectPbc3d::VecRefType& rhs) const
    {
        float x = dist1d(0, lhs[0], rhs[0]);
        float y = dist1d(1, lhs[1], rhs[1]);
        float z = dist1d(2, lhs[2], rhs[2]);

        return x*x + y*y + z*z;
    }


} // namespace libmd
