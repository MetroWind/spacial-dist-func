#ifndef SDF_CONFIG_H
#define SDF_CONFIG_H

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>

#include "trojectory.h"

namespace sdf
{
    struct Parameters
    {
        libmd::AtomIdentifier Anchor;
        libmd::AtomIdentifier AtomX;
        libmd::AtomIdentifier AtomXY;
        float Distance;
        float SliceThickness;
        std::unordered_set<libmd::AtomIdentifier> OtherAtoms;
    };

    struct RuntimeConfig
    {
        static RuntimeConfig read(std::istream& s);
        static RuntimeConfig readFromFile(const std::string& s);

        std::string XtcFile;
        std::string GroFile;
        std::vector<Parameters> Params;
    };

}

#endif
