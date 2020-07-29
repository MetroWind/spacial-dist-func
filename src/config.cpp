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

#include <cstdlib>
#include <fstream>

#include "config.h"

namespace sdf
{
    using namespace libmd;

    HCenter :: HCenter(const std::string& which)
    {
        if(which == "anchor")
        {
            Center = ANCHOR;
        }
        else if(which == "x")
        {
            Center = X;
        }
        else if(which == "xy")
        {
            Center = XY;
        }
        else
        {
            throw std::invalid_argument(std::string("Unknown H center type: ") + which);
        }
    }

    RuntimeConfig RuntimeConfig :: read(std::istream& s)
    {
        RuntimeConfig Config;

        std::string Buffer;

        // XTC filename
        std::getline(s, Buffer);
        Config.XtcFile = Buffer;
        // GRO filename
        std::getline(s, Buffer);
        Config.GroFile = Buffer;
        // Skip first +++
        std::getline(s, Buffer);

        do
        {
            Parameters CurrentParams;
            // Anchor Atom
            std::getline(s, Buffer);
            CurrentParams.Anchor = Buffer;
            // X Atom
            std::getline(s, Buffer);
            CurrentParams.AtomX = Buffer;
            // XY Atom
            std::getline(s, Buffer);
            CurrentParams.AtomXY = Buffer;
            // Distance
            std::getline(s, Buffer);
            CurrentParams.Distance = std::atof(Buffer.c_str());
            // Thickness
            std::getline(s, Buffer);
            CurrentParams.SliceThickness = std::atof(Buffer.c_str());

            // The config file has the ability to explicitly specify
            // the atoms that should be included in the calculation.
            // But I later learnt that this ability is not needed. All
            // atoms are included except for the atoms that share name
            // with the anchor, x atom, and xy atom.
            //
            // I still retain the code to read them here, but the
            // values are not used.
            std::getline(s, Buffer);
            s.peek();
            while(Buffer != "+++" && !s.eof())
            {
                CurrentParams.OtherAtoms.insert(strip(Buffer));
                std::getline(s, Buffer);
                s.peek();
            }
            if(Buffer != "+++")
            {
                CurrentParams.OtherAtoms.insert(strip(Buffer));
            }

            Config.Params.emplace_back(std::move(CurrentParams));
        } while(!s.eof());

        return Config;
    }


    RuntimeConfig RuntimeConfig :: readFromFile(const std::string& s)
    {
        std::ifstream f(s.c_str());
        auto p = read(f);
        f.close();
        return p;
    }
}
