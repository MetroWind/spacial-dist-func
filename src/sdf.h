// -*- mode: c++; -*-
#ifndef SDF_SDF_H
#define SDF_SDF_H

#include <type_traits>
#include <unordered_set>

#include "utils.h"
#include "trojectory.h"
#include "pbc.h"
#include "config.h"

namespace libmd
{
    template<class FrameType>
    void wrapFrame(const AtomIdentifier& anchor_name, FrameType& frame)
    {
        static_assert(std::is_same<FrameType, Trojectory>::value ||
                      std::is_same<FrameType, TrojectorySnapshot>::value,
                      "FrameType can only be either Trojectory or "
                      "TrojectorySnapshot");

        const auto& Anchor = frame.vec(anchor_name);
        const auto& BoxDim = frame.meta().BoxDim;
        const RectPbc3d Pbc(BoxDim[0][0], BoxDim[1][1], BoxDim[2][2]);

        for(int32_t i = 0; i < frame.meta().AtomCount; i++)
        {
            Pbc.wrapVec(Anchor, frame.vec(i));
        }
    }

    template<class VecType, class FrameType>
    void shiftFrame(const VecType& by, FrameType& frame)
    {
        static_assert(std::is_same<FrameType, Trojectory>::value ||
                      std::is_same<FrameType, TrojectorySnapshot>::value,
                      "FrameType can only be either Trojectory or "
                      "TrojectorySnapshot");
        for(int32_t i = 0; i < frame.meta().AtomCount; i++)
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
        for(int32_t i = 0; i < frame.meta().AtomCount; i++)
        {
            frame.vec(i) = rot * frame.vec(i);
        }
    }

    // Return a rotation matrix, which would rotate vec to +x, and
    // in_xy to somewhere in the xy plane; in_xy x vec would point to
    // -z (which means the y component of in_xy > 0).
    Eigen::Matrix3f rotateToAlignX(const VecRefType& vec,
                                   const VecRefType& to_xy);
} // namespace libmd

namespace sdf
{
    libmd::TrojectorySnapshot prepareFrame(
        const Parameters& params, const libmd::Trojectory& frame);

    void run(const RuntimeConfig& config);

} // namespace sdf

#endif
