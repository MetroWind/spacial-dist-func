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

#include <iostream>
#include <string>
#include <thread>

#include <getopt.h>

#include "sdf.h"

void usage(const std::string& prog_name)
{
    std::cout << "Usage: " << prog_name << " [OPTIONS] INPUT\n";
}

void help(const std::string& prog_name)
{
    usage(prog_name);
    std::cout << "\nOptions:\n"
        "\n"
"-h, --help                     Show this messages\n\n"
"-t N, --threads N              Number of threads to use. Default:\n"
"                               Number of CPU cores.\n\n"
"-d X, --distance-cutoff X      The cutoff distance when selecting\n"
"                               nearest atoms, in XTC native unit.\n"
"                               Override the value specified in the\n"
"                               input file.\n\n"
"-s X, --slice-thickness X      The thickness of slice of projection,\n"
"                               in XTC native unit. Override the value\n"
"                               specified in the input file.\n\n"
"-r N, --resolution X           The resolution of the histogram grid.\n"
"                               This is the number of grid cells in\n"
"                               each direction. Default: 40\n\n"
"--hist-range X                 The size of the 2D historgram box,\n"
"                               denoted by a ratio of the box size in\n"
"                               the XTC file (only the x and y size\n"
"                               matter). For example if this is set to\n"
"                               0.1, and the XTC box size is 3x4, the\n"
"                               histogram box is then 0.3x0.4. This\n"
"                               box is divided into resolution^2\n"
"                               number of grid cells. Default: 0.1\n\n"
"-p, --progress                 Show a “progress bar”.\n\n"
        ;
}

int main(int argc, char** argv)
{
    std::string ProgName(argv[0]);
    sdf::Parameters Params;
    size_t Threads = 0;
    float Distance = 0.0;
    float Thickness = 0.0;
    size_t Resolution = 40;
    float HistRange = 0.1;
    bool Progress = false;

    {
        static struct option Options[] = {
            { "help", no_argument, nullptr, 'h' },
            { "threads", required_argument, nullptr, 't' },
            { "distance-cutoff", required_argument, nullptr, 'd' },
            { "slice-thickness", required_argument, nullptr, 's' },
            { "resolution", required_argument, nullptr, 'r' },
            { "hist-range", required_argument, nullptr, 'H' },
            { "progress", required_argument, nullptr, 'p' },
            { nullptr, 0, nullptr, 0 }
        };

        int ch;
        while ((ch = getopt_long(argc, argv, "ht:d:s:r:p", Options, nullptr)) != -1)
        {
            switch (ch)
            {
            case 'h':
                help(ProgName);
                return 0;
            case 't':
                Threads = std::atoi(optarg);
                break;
            case 'd':
                Distance = std::atof(optarg);
                break;
            case 's':
                Thickness = std::atof(optarg);
                break;
            case 'r':
                Resolution = std::atoi(optarg);
                break;
            case 'H':
                HistRange = std::atof(optarg);
                break;
            case 'p':
                Progress = true;
                break;
            default:
                usage(ProgName);
                return -1;
            }
        }
        argc -= optind;
        argv += optind;
    }

    if(argc != 1)
    {
        usage(ProgName);
        return -1;
    }

    std::string InputFile = argv[0];

    auto Config = sdf::RuntimeConfig::readFromFile(InputFile);
    if(Threads == 0)
    {
        Config.ThreadCount = std::thread::hardware_concurrency();
    }
    else
    {
        Config.ThreadCount = Threads;
    }

    if(Distance > 0.0)
    {
        for(auto& Params: Config.Params)
        {
            Params.Distance = Distance;
        }
    }
    if(Thickness > 0.0)
    {
        for(auto& Params: Config.Params)
        {
            Params.SliceThickness = Thickness;
        }
    }
    Config.Resolution = Resolution;
    if(HistRange > 0.0)
    {
        Config.HistRange = HistRange;
    }
    Config.Progress = Progress;

    sdf::run(Config);

    return 0;
}
