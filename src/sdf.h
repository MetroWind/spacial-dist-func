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
#include <atomic>
#include <thread>
#include <mutex>
#include <iomanip>
#include <stdexcept>

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
    struct DistCountTraits
    {
        using ValueType = uint64_t;
        static const ValueType Zero;
    };

    struct DistChargeTraits
    {
        using ValueType = int64_t;
        static const ValueType Zero;
    };

    struct DistDetailedCountTraits
    {
        using ValueType = std::unordered_map<std::string, uint64_t>;
        static const ValueType Zero;
    };

    template <class DistTraits>
    class Distribution2
    {
    public:
        inline void buildGrid();
        inline void cornerLow(float x, float y);
        inline void cornerHigh(float x, float y);
        inline void resolution(size_t n);
        static inline typename DistTraits::ValueType
        deltaFromAtom(const libmd::AtomIdentifier& atom, const RuntimeConfig& config);
        inline void delta(float x, float y, typename DistTraits::ValueType d);
        inline void add(float x, float y);
        inline typename DistTraits::ValueType value(size_t ix, size_t iy) const;

        inline std::string prettyPrint() const;
        inline std::string jsonMesh(bool avg_over_frames) const;
        inline std::string jsonValue(
            const typename DistTraits::ValueType& x, bool avg_over_frames) const;

        inline void addSpecial(const std::string& name, std::array<float, 2> coord);

        size_t FrameCount;

    private:
        inline size_t index(size_t ix, size_t iy) const;
        inline size_t index(float x, float y) const;

        std::vector<typename DistTraits::ValueType> Count;
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
        const auto& AtomXY = frame.vec(params.AtomXY);
        const auto& Anchor = frame.vec(params.Anchor);

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

                switch(params.Center.type())
                {
                case HCenter::X:
                    return Pbc.dist(AtomX, vec) < params.Distance;
                case HCenter::XY:
                    return Pbc.dist(AtomXY, vec) < params.Distance;
                case HCenter::ANCHOR:
                    return Pbc.dist(Anchor, vec) < params.Distance;
                default:
                    throw std::runtime_error("Invalid center");
                }

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

    // A Dummy type to force the use of specialized deltaFromAtom();
    template<typename T> struct ForceSpecialization: public std::false_type {};

    template <class DistTraits>
    inline typename DistTraits::ValueType Distribution2<DistTraits> ::
    deltaFromAtom(const libmd::AtomIdentifier& atom, const RuntimeConfig& config)
    {
        UNUSED(atom); UNUSED(config);
        static_assert(ForceSpecialization<DistTraits>::value,
                      "Need implementation of deltaFromAtom()");
    }

    template <>
    inline typename DistCountTraits::ValueType Distribution2<DistCountTraits> ::
    deltaFromAtom(const libmd::AtomIdentifier& atom, const RuntimeConfig& config)
    {
        UNUSED(atom); UNUSED(config);
        return 1;
    }

    template <>
    inline typename DistChargeTraits::ValueType Distribution2<DistChargeTraits> ::
    deltaFromAtom(const libmd::AtomIdentifier& atom, const RuntimeConfig& config)
    {
        return config.AtomProperties.at(atom.Name).Charge;
    }

    template <>
    inline typename DistDetailedCountTraits::ValueType Distribution2<DistDetailedCountTraits> ::
    deltaFromAtom(const libmd::AtomIdentifier& atom, const RuntimeConfig& _)
    {
        UNUSED(_);
        typename DistDetailedCountTraits::ValueType Delta;
        Delta[atom.Name] = 1;
        return Delta;
    }

    template <class DistTraits>
    inline Distribution2<DistTraits> run(const RuntimeConfig& config)
    {
        libmd::Trajectory t;
        t.open(config.XtcFile, config.GroFile);

        Distribution2<DistTraits> Result;
        if(config.AbsoluteHistRange)
        {
            Result.cornerLow(-0.5 * config.HistRange,
                             -0.5 * config.HistRange);
            Result.cornerHigh(0.5 * config.HistRange,
                              0.5 * config.HistRange);
        }
        else
        {
            Result.cornerLow(-(t.meta().BoxDim[0][0]) * config.HistRange,
                             -(t.meta().BoxDim[1][1]) * config.HistRange);
            Result.cornerHigh(t.meta().BoxDim[0][0] * config.HistRange,
                              t.meta().BoxDim[1][1] * config.HistRange);
        }
        Result.resolution(config.Resolution);
        Result.buildGrid();

        t.nextFrame();
        {
            auto FirstFrame = prepareFrame(config.Params[0], t);
            auto AnchorName = config.Params[0].Anchor.toStr();
            auto Anchor = FirstFrame.extraAtoms().at(config.Params[0].Anchor);
            auto AtomXName = config.Params[0].AtomX.toStr();
            auto AtomX = FirstFrame.extraAtoms().at(config.Params[0].AtomX);
            auto AtomXYName = config.Params[0].AtomXY.toStr();
            auto AtomXY = FirstFrame.extraAtoms().at(config.Params[0].AtomXY);

            Result.addSpecial(AnchorName, {Anchor[0], Anchor[1]});
            Result.addSpecial(AtomXName, {AtomX[0], AtomX[1]});
            Result.addSpecial(AtomXYName, {AtomXY[0], AtomXY[1]});
        }

        // Make sure the atoms specified in the input exist.
        for(const auto& param: config.Params)
        {
            if(!t.hasAtom(param.AtomX))
            {
                throw std::runtime_error(std::string("Unknown atom: ") +
                                         param.AtomX.toStr());
            }
            if(!t.hasAtom(param.AtomXY))
            {
                throw std::runtime_error(std::string("Unknown atom: ") +
                                         param.AtomXY.toStr());
            }
            if(!t.hasAtom(param.Anchor))
            {
                throw std::runtime_error(std::string("Unknown atom: ") +
                                         param.Anchor.toStr());
            }
        }

        t.close();
        t.clear();
        t.open(config.XtcFile, config.GroFile);

        std::mutex HistLock;
        std::mutex FrameLock;
        std::vector<std::thread> Threads;

        for(size_t i = 0; i < config.ThreadCount; i++)
        {
            Threads.emplace_back(std::thread([&]()
            {
                while(true)
                {
                    FrameLock.lock();
                    if(!t.nextFrame())
                    {
                        FrameLock.unlock();
                        break;
                    }
                    const auto Frame = t.snapshot();
                    FrameLock.unlock();

                    for(const auto& Params: config.Params)
                    {
                        auto Frame = prepareFrame(Params, t);
                        HistLock.lock();
                        const auto Shot = Frame.snapshot();
                        for(size_t AtomIdx = 0; AtomIdx < Shot.size(); AtomIdx++)
                        {
                            const auto& Vec = Shot.vec(AtomIdx);
                            const auto& Atom = Shot.atomId(AtomIdx);
                            try
                            {
                                Result.delta(Vec[0], Vec[1],
                                             Distribution2<DistTraits>::
                                             deltaFromAtom(Atom, config));
                            }
                            catch(std::out_of_range)
                            {
                            }
                        }
                        HistLock.unlock();
                    }
                    if(config.Progress)
                    {
                        std::cerr << "." << std::flush;
                    }
                }
            }));
        }

        for(auto& Thread: Threads)
        {
            Thread.join();
        }

        if(config.Progress) { std::cerr << std::endl; }
        t.close();
        Result.FrameCount = t.countFrames();
        return Result;
    }

    template <class DistTraits>
    inline void Distribution2<DistTraits> :: cornerLow(float x, float y)
    {
        CornerLow[0] = x;
        CornerLow[1] = y;
    }

    template <class DistTraits>
    inline void Distribution2<DistTraits> :: cornerHigh(float x, float y)
    {
        CornerHigh[0] = x;
        CornerHigh[1] = y;
    }

    template <class DistTraits>
    inline void Distribution2<DistTraits> :: resolution(size_t n)
    {
        Resolution = n;
    }

    template <class DistTraits>
    inline void Distribution2<DistTraits> :: buildGrid()
    {
        Count.clear();
        Count.resize(Resolution * Resolution, DistTraits::Zero);
    }

    template <class DistTraits>
    inline void Distribution2<DistTraits> :: delta(
        float x, float y, typename DistTraits::ValueType d)
    {
        Count.at(index(x, y)) += d;
    }

    template <>
    inline void Distribution2<DistDetailedCountTraits> :: delta(
        float x, float y, typename DistDetailedCountTraits::ValueType d)
    {
        const auto Pair = std::begin(d);
        auto& Map = Count.at(index(x, y));
        auto Existing = Map.find(Pair->first);
        if(Existing == std::end(Map))
        {
            Map[Pair -> first] = Pair -> second;
        }
        else
        {
            Map[Pair -> first] += Pair -> second;
        }
    }

    template <class DistTraits>
    inline void Distribution2<DistTraits> :: add(float x, float y)
    {
        delta(x, y, 1);
    }

    template <class DistTraits>
    inline typename DistTraits::ValueType Distribution2<DistTraits> ::
    value(size_t ix, size_t iy) const
    {
        return Count.at(index(ix, iy));
    }

    template <class DistTraits>
    inline size_t Distribution2<DistTraits> :: index(size_t ix, size_t iy) const
    {
        if(ix >= Resolution || iy >= Resolution)
        {
            throw std::out_of_range("Grid index out of range");
        }
        return ix * Resolution + iy;
    }

    template <class DistTraits>
    inline size_t Distribution2<DistTraits> :: index(float x, float y) const
    {
        if(x < CornerLow[0] || x >= CornerHigh[0] || y < CornerLow[1] ||
           y >= CornerHigh[1])
        {
            throw std::out_of_range("Grid coordinate out of range");
        }

        size_t ix = (x - CornerLow[0]) / (CornerHigh[0] - CornerLow[0]) *
            float(Resolution);
        size_t iy = (y - CornerLow[1]) / (CornerHigh[1] - CornerLow[1]) *
            float(Resolution);
        return index(ix, iy);
    }

    template <class DistTraits>
    inline void Distribution2<DistTraits> :: addSpecial(
        const std::string& name, std::array<float, 2> coord)
    {
        Specials[name] = coord;
    }

    template <class DistTraits>
    inline std::string Distribution2<DistTraits> :: prettyPrint() const
    {
        std::stringstream Formatter;
        for(size_t y = 0; y < Resolution; y++)
        {
            for(size_t x = 0; x < Resolution; x++)
            {
                Formatter << std::setw(8) << value(x, y);
            }
            Formatter << "\n";
        }
        return Formatter.str();
    }

    template <class DistTraits>
    inline std::string Distribution2<DistTraits> :: jsonValue(
        const typename DistTraits::ValueType& x, bool avg_over_frames) const
    {
        std::stringstream Formatter;
        if(avg_over_frames)
        {
            Formatter << float(x) / float(FrameCount);
        }
        else
        {
            Formatter << x;
        }
        return Formatter.str();
    }

    template <>
    inline std::string Distribution2<DistDetailedCountTraits> :: jsonValue(
        const typename DistDetailedCountTraits::ValueType& x, bool avg_over_frames) const
    {
        if(x.empty())
        {
            return "{}";
        }
        std::stringstream Formatter;
        Formatter << "{\n";
        for(const auto& Pair: x)
        {
            Formatter << "  \"" << Pair.first << "\": ";
            if(avg_over_frames)
            {
                Formatter << float(Pair.second) / float(FrameCount);
            }
            else
            {
                Formatter << float(Pair.second);
            }
            Formatter << ",\n";
        }
        Formatter.seekp(-2, Formatter.cur);
        Formatter << "} ";
        return Formatter.str();
    }

    template <class DistTraits>
    inline std::string Distribution2<DistTraits> :: jsonMesh(bool avg_over_frames) const
    {
        std::stringstream Formatter;
        Formatter << "{\n";
        Formatter << "\"c\": [";
        for(size_t y = 0; y < Resolution; y++)
        {
            Formatter << "[";
            for(size_t x = 0; x < Resolution; x++)
            {
                Formatter << jsonValue(value(x, y), avg_over_frames);
                if(x < Resolution - 1)
                {
                    Formatter << ", ";
                }
            }
            Formatter << "]";
            if(y < Resolution - 1)
            {
                Formatter << ",";
                Formatter << "\n";
            }
        }
        Formatter << "],\n";

        // Output the mesh for the grid. The structure of the mesh is
        // documented at
        // https://matplotlib.org/api/_as_gen/matplotlib.axes.Axes.pcolormesh.html#matplotlib.axes.Axes.pcolormesh
        Formatter << "\"x\": [";
        float XSize = CornerHigh[0] - CornerLow[0];
        float CellSizeX = XSize / float(Resolution);
        for(size_t iy = 0; iy <= Resolution; iy++)
        {
            Formatter << "[";
            for(size_t ix = 0; ix <= Resolution; ix++)
            {
                Formatter << float(ix) * CellSizeX + CornerLow[0];
                if(ix < Resolution)
                {
                    Formatter << ", ";
                }
            }
            Formatter << "]";
            if(iy < Resolution)
            {
                Formatter << ",";
                Formatter << "\n";
            }
        }
        Formatter << "],\n";

        Formatter << "\"y\": [";
        float YSize = CornerHigh[1] - CornerLow[1];
        float CellSizeY = YSize / float(Resolution);
        for(size_t iy = 0; iy <= Resolution; iy++)
        {
            Formatter << "[";
            float y = float(iy) * CellSizeY + CornerLow[1];
            for(size_t ix = 0; ix <= Resolution; ix++)
            {
                Formatter << y;
                if(ix < Resolution)
                {
                    Formatter << ", ";
                }
            }
            Formatter << "]";
            if(iy < Resolution)
            {
                Formatter << ",";
                Formatter << "\n";
            }
        }
        Formatter << "],\n";

        Formatter << "\"specials\": {";
        for(const auto& SpecialPair: Specials)
        {
            Formatter << "\"" << SpecialPair.first << "\": [";
            Formatter << SpecialPair.second[0] << ", "
                      << SpecialPair.second[1] << "],\n";
        }
        Formatter.seekp(-2, Formatter.cur);
        Formatter << "\n}}\n";

        return Formatter.str();
    }

} // namespace sdf

#endif
