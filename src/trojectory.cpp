#include <fstream>
#include <array>
#include <sstream>
#include <cstring>

#include "trojectory.h"

namespace libmd
{
    namespace
    {
        // Equivalent of Python’s str.strip().
        std::string strip(const std::string& str,
                          const std::string& whitespace = " \t\n")
        {
            const auto StrBegin = str.find_first_not_of(whitespace);
            if (StrBegin == std::string::npos)
                return ""; // no content

            const auto StrEnd = str.find_last_not_of(whitespace);
            const auto StrRange = StrEnd - StrBegin + 1;

            return str.substr(StrBegin, StrRange);
        }

        // Extract names from .gro file. Format spec:
        // http://manual.gromacs.org/current/reference-manual/file-formats.html#gro.
        std::vector<std::string> extractAtomNames(std::ifstream& file)
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

            std::vector<std::string> Result;
            for(size_t i = 0; i < Count; i++)
            {
                file.getline(Buffer.data(), Buffer.size());
                // The name starts at column 11, and is 5 chars long.
                std::string Name(Buffer.data() + 10, Buffer.data() + 15);
                Result.push_back(strip(Name));
            }

            return Result;
        }

    } // namespace

    void Trojectory :: open(const std::string& xtc_path,
                            const std::string& gro_path)
    {
        f.open(xtc_path.c_str());

        std::ifstream GroFile(gro_path);
        AtomNames = extractAtomNames(GroFile);

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

    bool Trojectory :: nextFrame()
    {
        if(f.eof())
        {
            return false;
        }

        f.readFrame(Data.data());
        return true;
    }

    void Trojectory :: close()
    {
        if(f.isOpen()) { f.close(); }
        Vecs.clear();
        Data.clear();
    }

} // namespace libmd
