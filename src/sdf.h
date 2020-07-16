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

#ifndef SDF_SDF_H
#define SDF_SDF_H

#include <type_traits>
#include <unordered_set>

#include "utils.h"
#include "trajectory.h"
#include "pbc.h"
#include "config.h"

namespace libmd
{
    template<class FrameType>
    void wrapFrame(const AtomIdentifier& anchor_name, FrameType& frame)
    {
        static_assert(std::is_same<FrameType, Trajectory>::value ||
                      std::is_same<FrameType, TrajectorySnapshot>::value,
                      "FrameType can only be either Trajectory or "
                      "TrajectorySnapshot");

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
        static_assert(std::is_same<FrameType, Trajectory>::value ||
                      std::is_same<FrameType, TrajectorySnapshot>::value,
                      "FrameType can only be either Trajectory or "
                      "TrajectorySnapshot");
        for(int32_t i = 0; i < frame.meta().AtomCount; i++)
        {
            frame.vec(i) += by;
        }
    }

    template<class FrameType>
    void rotateFrame(const Eigen::Ref<Eigen::Matrix3f>& rot, FrameType& frame)
    {
        static_assert(std::is_same<FrameType, Trajectory>::value ||
                      std::is_same<FrameType, TrajectorySnapshot>::value,
                      "FrameType can only be either Trajectory or "
                      "TrajectorySnapshot");
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
    class Distribution2
    {
    public:
        void buildGrid();
        void cornerLow(float x, float y);
        void cornerHigh(float x, float y);
        void resolution(size_t n);
        void add(float x, float y);
        uint32_t count(size_t ix, size_t iy) const;

        std::string prettyPrint() const;
        std::string jsonMesh(bool avg_over_frames) const;

        void addSpecial(const std::string& name, std::array<float, 2> coord);

        size_t FrameCount;

    private:
        size_t index(size_t ix, size_t iy) const;
        size_t index(float x, float y) const;

        std::vector<uint32_t> Count;
        std::array<float, 2> CornerLow;
        std::array<float, 2> CornerHigh;
        size_t Resolution;    // Number of grid cells in each direction
        std::unordered_map<std::string, std::array<float, 2>> Specials;
    };

    class PreparedFrame
    {
    public:
        PreparedFrame() = delete;
        PreparedFrame(const libmd::TrajectorySnapshot& frame)
                : Frame(frame) {}
        void addExtra(const libmd::AtomIdentifier& id,
                      const libmd::VecRefType& vec)
        {
            ExtraAtoms[id] = vec;
        }

        const std::unordered_map<libmd::AtomIdentifier, Eigen::Vector3f>
        extraAtoms() const
        {
            return ExtraAtoms;
        }

        const libmd::TrajectorySnapshot& snapshot() const
        {
            return Frame;
        }

    private:
    std::unordered_map<libmd::AtomIdentifier, Eigen::Vector3f> ExtraAtoms;
    const libmd::TrajectorySnapshot Frame;
    };

    template <class FrameType>
    PreparedFrame prepareFrame(
        const Parameters& params, const FrameType& frame)
    {
        static_assert(std::is_same<FrameType, libmd::Trajectory>::value ||
                      std::is_same<FrameType, libmd::TrajectorySnapshot>::value,
                      "FrameType can only be either Trajectory or "
                      "TrajectorySnapshot");

        const auto& AtomX = frame.vec(params.AtomX);
        const auto& BoxDim = frame.meta().BoxDim;
        const libmd::RectPbc3d Pbc(BoxDim[0][0], BoxDim[1][1], BoxDim[2][2]);

        // Filter by distance
        auto Frame = filterFrame(
            frame,
            [&](const libmd::AtomIdentifier& name, size_t _, const auto& vec)
            {
                UNUSED(_);
                if(name == params.Anchor || name == params.AtomX ||
                   name == params.AtomXY)
                {
                    return true;
                }

                if(name.Res == params.Anchor.Res)
                {
                    return false;
                }

                // if(params.OtherAtoms.find(name) == std::end(params.OtherAtoms))
                // {
                //     return false;
                // }

                return Pbc.dist(AtomX, vec) < params.Distance;
            });

        // std::cout << Frame.debugString() << std::endl;

        wrapFrame(params.AtomX, Frame);
        Eigen::Vector3f ShiftBy = -(Frame.vec(params.Anchor));
        shiftFrame(ShiftBy, Frame);
        auto Rot = libmd::rotateToAlignX(Frame.vec(params.AtomX),
                                  Frame.vec(params.AtomXY));
        rotateFrame(Rot, Frame);
        // std::cerr << Frame.debugString() << std::endl;

        // Ditch the anchors, x atoms, and xy atoms, and take a slice
        // at the XY plane.
        float HalfThickness = params.SliceThickness * 0.5;
        auto Snap = filterFrame(
            Frame,
            [&](const libmd::AtomIdentifier& id, auto _1, const auto& pos)
            {
                UNUSED(_1);
                return !(id.Name == params.Anchor.Name ||
                         id.Name == params.AtomX.Name ||
                         id.Name == params.AtomXY.Name) &&
                    std::fabs(pos[2]) <= HalfThickness;
            });
        PreparedFrame Result(std::move(Snap));
        Result.addExtra(params.Anchor, Frame.vec(params.Anchor));
        Result.addExtra(params.AtomX, Frame.vec(params.AtomX));
        Result.addExtra(params.AtomXY, Frame.vec(params.AtomXY));
        // for(size_t ih = 1; ih <= 14; ih++)
        // {
        //     std::stringstream Formatter;
        //     Formatter << "1+H" << ih;
        //     const auto Name = Formatter.str();
        //     if(Frame.hasAtom(Name))
        //     {
        //         Result.addExtra(Name, Frame.vec(Name));
        //     }
        // }

        return Result;
    }

    Distribution2 run(const RuntimeConfig& config);

} // namespace sdf

#endif
