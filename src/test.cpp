#include <iostream>
#include <string>

#include <getopt.h>

#include "sdf.h"

void usage(const std::string& prog_name)
{
    std::cout << prog_name << " [OPTIONS] INPUT\n";
}

void help(const std::string& prog_name)
{
    usage(prog_name);
    std::cout << "Options:\n"
"\n"
"-h, --help                     Show this messages\n"
        ;
}

int main(int argc, char** argv)
{
    std::string ProgName(argv[0]);
    sdf::Parameters Params;

    {
        static struct option Options[] = {
            { "help", no_argument, nullptr, 'h' },
            { nullptr, 0, nullptr, 0 }
        };

        int ch;
        while ((ch = getopt_long(argc, argv, "h", Options, nullptr)) != -1) {
            switch (ch) {
            case 'h':
                help(ProgName);
                return 0;
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
    sdf::run(Config);

    return 0;
}
