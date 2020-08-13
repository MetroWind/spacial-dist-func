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

#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include <getopt.h>


#include "sdf.h"


void handler(int sig)
{
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    std::cerr << "Error: signal " << sig << std::endl;
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(128 + sig);
}


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
"    Number of CPU cores.\n\n"
"-d X, --distance-cutoff X      The cutoff distance when selecting\n"
"    nearest atoms, in XTC native unit. Override the value specified in\n"
"    the input file.\n\n"
"-s X, --slice-thickness X      The thickness of slice of projection,\n"
"    in XTC native unit. Override the value specified in the input\n"
"    file.\n\n"
"-r N, --resolution N           The resolution of the histogram grid.\n"
"    This is the number of grid cells in each direction. Default: 40\n\n"
"--hist-range X                 The size of the 2D historgram region,\n"
"    denoted by a ratio of the box size in the XTC file (only the x and\n"
"    y size matter). For example if this is set to 0.1, and the XTC box\n"
"    size is 3x4, the histogram box is then 0.3x0.4. This box is\n"
"    divided into resolution^2 number of grid cells. Default: 0.1\n\n"
"--hist-range-abs X             The absolute size of the 2D historgram\n"
"    region. This is mutually exclusive with --hist-range. Default: use\n"
"    --hist-range 0.1\n\n"
"-p, --progress                 Show a “progress bar”.\n\n"
"-a, --average                  Average the result over number of\n"
"    frames.\n\n"
"--center TYPE                  The type of center atom to calculate\n"
"    distance from. Valid types are 'anchor', 'x', and 'xy'. Default:\n"
"    x.\n\n"
"--measure TYPE                 The quantity of which to make\n"
"    distribution. Valid arguments are 'count', 'charge', and\n"
"    'count-per-atom'. Default: count.\n\n"
        ;
}

int main(int argc, char** argv)
{
    signal(SIGSEGV, handler);
    signal(SIGABRT, handler);

    std::string ProgName(argv[0]);
    sdf::Parameters Params;
    size_t Threads = 0;
    float Distance = 0.0;
    float Thickness = 0.0;
    size_t Resolution = 40;
    float HistRange = 0.1;
    bool AbsoluteHistRange = false;
    bool Progress = false;
    bool Average = false;
    sdf::HCenter CenterType("x");
    int MeasureSpecified = 0;
    std::string Measure("count");
    const std::unordered_set<std::string> ValidMeasures =
        {"count", "charge", "count-per-atom"};

    {
        static struct option Options[] = {
            { "help", no_argument, nullptr, 'h' },
            { "threads", required_argument, nullptr, 't' },
            { "distance-cutoff", required_argument, nullptr, 'd' },
            { "slice-thickness", required_argument, nullptr, 's' },
            { "resolution", required_argument, nullptr, 'r' },
            { "hist-range", required_argument, nullptr, 'n' },
            { "hist-range-abs", required_argument, nullptr, 'N' },
            { "progress", no_argument, nullptr, 'p' },
            { "average", no_argument, nullptr, 'a' },
            { "center", required_argument, nullptr, 'c' },
            { "measure", required_argument, &MeasureSpecified, 1},
            { nullptr, 0, nullptr, 0 }
        };

        int ch;
        while ((ch = getopt_long(argc, argv, "ht:d:s:r:pa", Options, nullptr)) != -1)
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
            case 'n':
                HistRange = std::atof(optarg);
                AbsoluteHistRange = false;
                break;
            case 'N':
                HistRange = std::atof(optarg);
                AbsoluteHistRange = true;
                break;
            case 'p':
                Progress = true;
                break;
            case 'a':
                Average = true;
                break;
            case 'c':
                CenterType = std::string(optarg);
                break;
            case 0:
                if(MeasureSpecified == 1)
                {
                    Measure = optarg;
                }
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

    if(ValidMeasures.find(Measure) == std::end(ValidMeasures))
    {
        std::cerr << "Invalid measure: " << Measure << std::endl;
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

    for(auto& Params: Config.Params)
    {
        Params.Center = CenterType;
    }

    Config.Resolution = Resolution;
    if(HistRange > 0.0)
    {
        Config.HistRange = HistRange;
        Config.AbsoluteHistRange = AbsoluteHistRange;
    }
    Config.Progress = Progress;
    Config.AverageOverFrameCount = Average;

    if(Measure == "count")
    {
        const auto Result = sdf::run<sdf::DistCountTraits>(Config);
        std::cout << Result.jsonMesh(Config.AverageOverFrameCount);
    }
    else if(Measure == "charge")
    {
        const auto Result = sdf::run<sdf::DistChargeTraits>(Config);
        std::cout << Result.jsonMesh(Config.AverageOverFrameCount);
    }
    else if(Measure == "count-per-atom")
    {
        const auto Result = sdf::run<sdf::DistDetailedCountTraits>(Config);
        std::cout << Result.jsonMesh(Config.AverageOverFrameCount);
    }
    else
    {
        std::cerr << "Shit happened!" << std::endl;
        return 1;
    }

    return 0;
}
