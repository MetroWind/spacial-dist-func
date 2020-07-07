// -*- mode: c++; -*-
#ifndef SDF_SDF_H
#define SDF_SDF_H

#include <type_traits>

#include "utils.h"
#include "trojectory.h"
#include "pbc.h"

namespace libmd
{
    template<class FrameType>
    TrojectorySnapshot findNearest(
        const std::string& name, const float cutoff, const FrameType& frame)
    {
        static_assert(std::is_same<FrameType, Trojectory>::value ||
                      std::is_same<FrameType, TrojectorySnapshot>::value,
                      "FrameType can only be either Trojectory or "
                      "TrojectorySnapshot");

        const auto& BoxDim = frame.meta().BoxDim;
        const RectPbc3d Pbc(BoxDim[0][0], BoxDim[1][1], BoxDim[2][2]);

        const auto& Anchor = frame.vec(name);

        auto Filter = [&](auto& _1, auto& _2, const V3Map& vec)
        {
            UNUSED(_1); UNUSED(_2);
            return Pbc.dist(vec, Anchor) < cutoff;
        };

        return frame.filter(Filter);
    }

    template<class FrameType>
    void wrapFrame(const std::string& anchor_name, FrameType& frame)
    {
        static_assert(std::is_same<FrameType, Trojectory>::value ||
                      std::is_same<FrameType, TrojectorySnapshot>::value,
                      "FrameType can only be either Trojectory or "
                      "TrojectorySnapshot");

        const auto& Anchor = frame.vec(anchor_name);
        const auto& BoxDim = frame.meta().BoxDim;
        const RectPbc3d Pbc(BoxDim[0][0], BoxDim[1][1], BoxDim[2][2]);

        for(size_t i = 0; i < frame.meta().AtomCount; i++)
        {
            Pbc.wrapVec(Anchor, frame.vec(i));
        }
    }

    template<class FrameType>
    void shiftFrame(const VecRefType& by, FrameType& frame)
    {
        static_assert(std::is_same<FrameType, Trojectory>::value ||
                      std::is_same<FrameType, TrojectorySnapshot>::value,
                      "FrameType can only be either Trojectory or "
                      "TrojectorySnapshot");
        for(size_t i = 0; i < frame.meta().AtomCount; i++)
        {
            frame.vec(i) += by;
        }
    }

    template<class FrameType>
    void rotateFrame(const Eigen::Ref<Eigen::Matrix3f>& rot, FrameType& frame)
    {
        static_assert(std::is_same<FrameType, Trojectory>::value ||
                      std::is_same<FrameType, TrojectorySnapshot>::value,
                      "FrameType can only be either Trojectory or "
                      "TrojectorySnapshot");
        for(size_t i = 0; i < frame.meta().AtomCount; i++)
        {
            frame.vec(i) = rot * frame.vec(i);
        }
    }

    // Return a rotation matrix, which would rotate vec2x to +x, and
    // in_xy to somewhere in the xy plaine; in_xy x vec2x would point
    // to +z.
    Eigen::Matrix3f rotateToAlignX(const VecRefType& vec2x,
                                   const VecRefType& in_xy)
    {
    }

}

namespace sdf
{

} // namespace libmd

#endif
