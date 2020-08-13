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

#ifndef SDF_CONFIG_H
#define SDF_CONFIG_H

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include "trajectory.h"

namespace sdf
{
    class HCenter
    {
    public:
        using Type = enum {ANCHOR, X, XY};

        HCenter() = default;
        HCenter(const HCenter&) = default;
        HCenter(const std::string& which);
        HCenter(Type which) : Center(which) {}

        HCenter& operator=(const HCenter&) = default;

        Type type() const { return Center; }

    private:
        Type Center = X;
    };

    struct Parameters
    {
        libmd::AtomIdentifier Anchor;
        libmd::AtomIdentifier AtomX;
        libmd::AtomIdentifier AtomXY;
        float Distance;
        float SliceThickness;
        std::unordered_set<libmd::AtomIdentifier> OtherAtoms;
        HCenter Center;
    };

    struct AtomProperty
    {
        int8_t Charge;
    };

    using AtomPropertyMap = std::unordered_map<std::string, AtomProperty>;

    struct RuntimeConfig
    {
        static RuntimeConfig read(std::istream& s);
        static RuntimeConfig readFromFile(const std::string& s);

#ifndef STUPID_UBUNTU
        void printXml() const;
#endif

        std::string XtcFile;
        std::string GroFile;
        std::vector<Parameters> Params;
        size_t Resolution = 40;
        float HistRange = 0.1;
        bool AbsoluteHistRange = false;
        bool Progress = false;
        size_t ThreadCount = 0;
        bool AverageOverFrameCount = false;
        AtomPropertyMap AtomProperties;
    };

}

#endif
