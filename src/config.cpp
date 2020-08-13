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
#include <iostream>

#include "pugixml.hpp"

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

        pugi::xml_document Input;
        if(!Input.load(s))
        {
            throw std::runtime_error("Failed to parse input");
        }
        Config.XtcFile = Input.child("sdf-run").child("input")
            .child("trajectory").text().as_string();
        Config.GroFile = Input.child("sdf-run").child("input")
            .child("structure").text().as_string();

        for(const auto& AtomProp: Input.child("atom-properties").children("atom"))
        {
            const char* NameRaw = AtomProp.attribute("name").as_string();
            if(NameRaw == nullptr)
            {
                throw std::runtime_error("Atom property with empty name");
            }
            const std::string Name(NameRaw);
            const auto& ChargeNode = AtomProp.child("charge");
            if(ChargeNode.empty())
            {
                continue;
            }

            AtomProperty Props;
            Props.Charge = ChargeNode.text().as_int();
            Config.AtomProperties[Name] = std::move(Props);
        }

        for(const auto& Basis: Input.child("sdf-run").child("bases")
                .children("basis"))
        {
            Parameters Params;
            Params.Anchor = strip(Basis.child("anchor").text().as_string());
            Params.AtomX = strip(Basis.child("x").text().as_string());
            Params.AtomXY = strip(Basis.child("xy").text().as_string());
            Params.Distance = std::atof(
                strip(Basis.child("search-radius").text().as_string()).c_str());
            Params.SliceThickness = std::atof(
                strip(Basis.child("thickness").text().as_string()).c_str());
            for(const auto& Exclude: Basis.child("excludes").children("exclude"))
            {
                Params.OtherAtoms.insert(strip(Exclude.text().as_string()));
            }
            Config.Params.emplace_back(std::move(Params));
        }

        return Config;
    }

    RuntimeConfig RuntimeConfig :: readLegacy(std::istream& s)
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

    void RuntimeConfig :: printXml() const
    {
        std::cout << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
        std::cout << "<sdf-run>" << std::endl;
        std::cout << "<input>" << std::endl;
        std::cout << "<trajectory>" << XtcFile << "</trajectory>" << std::endl;
        std::cout << "<structure>" << GroFile << "</structure>" << std::endl;
        std::cout << "</input>" << std::endl;
        std::cout << "<config/>" << std::endl;
        std::cout << "<bases>" << std::endl;
        for(const auto& Param: Params)
        {
            std::cout << "<basis>" << std::endl;
            std::cout << "<anchor>" << Param.Anchor << "</anchor>" << std::endl;
            std::cout << "<x>" << Param.AtomX << "</x>" << std::endl;
            std::cout << "<xy>" << Param.AtomXY << "</xy>" << std::endl;
            std::cout << "<search-radius>" << Param.Distance << "</search-radius>" << std::endl;
            std::cout << "<thickness>" << Param.SliceThickness << "</thickness>" << std::endl;
            std::cout << "<excludes/>" << std::endl;
            std::cout << "</basis>" << std::endl;
        }
        std::cout << "</bases>\n</sdf-run>" << std::endl;
    }

}
