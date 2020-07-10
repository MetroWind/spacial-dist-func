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

#include <fstream>
#include <array>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "trajectory.h"

namespace libmd
{
    namespace
    {
        // Extract names from .gro file. Format spec:
        // http://manual.gromacs.org/current/reference-manual/file-formats.html#gro.
        std::vector<AtomIdentifier> extractAtomIds(std::ifstream& file)
        {
            std::array<char, 1024> Buffer;
            // Skip 1st line
            file.getline(Buffer.data(), Buffer.size());

            // 2nd line is number of atoms
            file.getline(Buffer.data(), Buffer.size());
            size_t Count = 0;
            {
                std::stringstream Reader(Buffer.data());
                Reader >> Count;
            }

            std::vector<AtomIdentifier> Result;
            for(size_t i = 0; i < Count; i++)
            {
                file.getline(Buffer.data(), Buffer.size());
                // The name starts at column 11, and is 5 chars long.
                const std::string Res(Buffer.data() + 0, Buffer.data() + 5);
                const std::string Name(Buffer.data() + 10, Buffer.data() + 15);

                AtomIdentifier Id(std::atoi(Res.c_str()), strip(Name));
                Result.emplace_back(std::move(Id));
            }

            return Result;
        }

    } // namespace

    AtomIdentifier :: AtomIdentifier(const std::string& s)
    {
        auto SepPos = s.find("+");
        Res = std::atoi(s.substr(0, SepPos).c_str());
        Name = s.substr(SepPos + 1, s.size() - SepPos - 1);
    }

    AtomIdentifier :: AtomIdentifier(const char s[])
    {
        size_t i;
        size_t SepPos = 0;
        for(i = 0; s[i] != '\0'; i++)
        {
            if(s[i] == '+')
            {
                Res = std::atoi(std::string(s, i).c_str());
                SepPos = i;
            }
        }
        Name = std::string(s + SepPos + 1, i - SepPos - 1);
    }

    TrajectorySnapshot :: TrajectorySnapshot(const TrajectorySnapshot& from)
            : AtomNames(from.AtomNames),
              AtomNamesReverse(from.AtomNamesReverse),
              Meta(from.Meta),
              Data(from.Data)
    {
        Vecs.clear();
        for(size_t i = 0; i < AtomNames.size(); i++)
        {
            Vecs.emplace_back(&(Data[i*3]));
        }
    }

    TrajectorySnapshot& TrajectorySnapshot :: operator=(TrajectorySnapshot from)
    {
        swap(*this, from);
        return *this;
    }

    std::string TrajectorySnapshot :: debugString() const
    {
        std::stringstream Formatter;
        for(size_t i = 0; i < static_cast<size_t>(meta().AtomCount); i++)
        {
            std::array<char, 128> Buffer;
            std::sprintf(Buffer.data(), "%5s %5.3f %5.3f %5.3f\n",
                         AtomNames[i].toStr().c_str(),
                         vec(i)[0], vec(i)[1], vec(i)[2]);
            Formatter << Buffer.data();
        }
        return Formatter.str();
    }

    void Trajectory :: open(const std::string& xtc_path,
                            const std::string& gro_path)
    {
        f.open(xtc_path.c_str());

        std::ifstream GroFile(gro_path);
        AtomNames = extractAtomIds(GroFile);

        Meta = f.readFrameMeta();
        if(static_cast<size_t>(Meta.AtomCount) != AtomNames.size())
        {
            throw std::runtime_error(
                "number of atoms does not align between XTC and GRO");
        }

        for(size_t i = 0; i < AtomNames.size(); i++)
        {
            AtomNamesReverse[AtomNames[i]] = i;
        }

        GroFile.close();
        Data.resize(Meta.AtomCount * 3); // 3D vector
        for(size_t i = 0; i < AtomNames.size(); i++)
        {
            Vecs.emplace_back(&(Data[i*3]));
        }
    }

    bool Trajectory :: nextFrame()
    {
        if(f.eof())
        {
            return false;
        }

        f.readFrame(Data.data());
        return true;
    }

    void Trajectory :: close()
    {
        if(f.isOpen()) { f.close(); }
    }

    std::string Trajectory :: debugString() const
    {
        std::stringstream Formatter;
        for(size_t i = 0; i < static_cast<size_t>(meta().AtomCount); i++)
        {
            std::array<char, 128> Buffer;
            std::sprintf(Buffer.data(), "%5s %5.3f %5.3f %5.3f\n",
                         AtomNames[i].toStr().c_str(),
                         vec(i)[0], vec(i)[1], vec(i)[2]);
            Formatter << Buffer.data();
        }
        return Formatter.str();
    }

} // namespace libmd
