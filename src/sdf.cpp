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

#include <atomic>
#include <thread>
#include <mutex>
#include <iomanip>
#include <stdexcept>

#include "sdf.h"
#include "utils.h"

using namespace libmd;

namespace libmd
{
    Eigen::Matrix3f rotateToAlignX(const VecRefType& vec,
                                   const VecRefType& to_xy)
    {
        Eigen::Matrix3f Rot1;
        {
            Eigen::Vector3f Axis;
            Axis[0] = 0.0f;
            Axis[1] = vec[2];
            Axis[2] = -vec[1];
            Axis /= Axis.norm();

            const float Cos = vec[0] / vec.norm();
            const float Sin = std::sqrt(1.0f - Cos*Cos);

            Rot1(0, 0) = Cos;
            Rot1(0, 1) = -Axis[2] * Sin;
            Rot1(0, 2) = Axis[1] * Sin;
            Rot1(1, 0) = Axis[2] * Sin;
            Rot1(1, 1) = Cos + Axis[1] * Axis[1] * (1.0f - Cos);
            Rot1(1, 2) = Axis[1] * Axis[2] * (1.0f - Cos);
            Rot1(2, 0) = -Axis[1] * Sin;
            Rot1(2, 1) = Axis[2] * Axis[1] * (1.0f - Cos);
            Rot1(2, 2) = Cos + Axis[2] * Axis[2] * (1.0f - Cos);
        }

        const Eigen::Vector3f Rotated = Rot1 * to_xy;
        Eigen::Matrix3f Rot2 = Eigen::Matrix3f::Zero();
        {
            float Cos = Rotated[1] / std::sqrt(Rotated[1] * Rotated[1] +
                                                     Rotated[2] * Rotated[2]);
            float Sin = std::sqrt(1.0f - Cos*Cos);
            if(Rotated[2] > 0)
            {
                Sin = -Sin;
            }

            Rot2(0, 0) = 1.0f;
            Rot2(1, 1) = Cos;
            Rot2(1, 2) = -Sin;
            Rot2(2, 1) = Sin;
            Rot2(2, 2) = Cos;
        }

        return Rot2 * Rot1;
    }
} // namespace libmd

namespace sdf
{
    using namespace libmd;

    Distribution2 run(const RuntimeConfig& config)
    {
        libmd::Trajectory t;
        t.open(config.XtcFile, config.GroFile);

        Distribution2 Result;
        Result.cornerLow(-(t.meta().BoxDim[0][0]) * config.HistRange,
                         -(t.meta().BoxDim[1][1]) * config.HistRange);
        Result.cornerHigh(t.meta().BoxDim[0][0] * config.HistRange,
                          t.meta().BoxDim[1][1] * config.HistRange);
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
                        for(const auto& Vec: Frame.snapshot().vecs())
                        {
                            try
                            {
                                Result.add(Vec[0], Vec[1]);
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
        std::cout << Result.jsonMesh(config.AverageOverFrameCount);
        Result.FrameCount = t.countFrames();
        return Result;

    }

    void Distribution2 :: cornerLow(float x, float y)
    {
        CornerLow[0] = x;
        CornerLow[1] = y;
    }

    void Distribution2 :: cornerHigh(float x, float y)
    {
        CornerHigh[0] = x;
        CornerHigh[1] = y;
    }

    void Distribution2 :: resolution(size_t n)
    {
        Resolution = n;
    }

    void Distribution2 :: buildGrid()
    {
        Count.clear();
        Count.resize(Resolution * Resolution, 0);
    }

    void Distribution2 :: add(float x, float y)
    {
        ++(Count.at(index(x, y)));
    }

    uint32_t Distribution2 :: count(size_t ix, size_t iy) const
    {
        return Count.at(index(ix, iy));
    }

    size_t Distribution2 :: index(size_t ix, size_t iy) const
    {
        if(ix >= Resolution || iy >= Resolution)
        {
            throw std::out_of_range("Grid index out of range");
        }
        return ix * Resolution + iy;
    }

    size_t Distribution2 :: index(float x, float y) const
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

    void Distribution2 :: addSpecial(
        const std::string& name, std::array<float, 2> coord)
    {
        Specials[name] = coord;
    }

    std::string Distribution2 :: prettyPrint() const
    {
        std::stringstream Formatter;
        for(size_t y = 0; y < Resolution; y++)
        {
            for(size_t x = 0; x < Resolution; x++)
            {
                Formatter << std::setw(8) << count(x, y);
            }
            Formatter << "\n";
        }
        return Formatter.str();
    }

    std::string Distribution2 :: jsonMesh(bool avg_over_frames) const
    {
        std::stringstream Formatter;
        Formatter << "{\n";
        Formatter << "\"c\": [";
        for(size_t y = 0; y < Resolution; y++)
        {
            Formatter << "[";
            for(size_t x = 0; x < Resolution; x++)
            {
                if(avg_over_frames)
                {
                    Formatter << float(count(x, y)) / float(FrameCount);
                }
                else
                {
                    Formatter << count(x, y);
                }

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
