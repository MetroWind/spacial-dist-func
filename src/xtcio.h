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

#ifndef SDF_XTCIO_H
#define SDF_XTCIO_H

#include <iostream>
#include <fstream>
#include <array>

#include "endian.h"

namespace libmd
{
    class XtcFile
    {
    public:
        static const int32_t MAGIC = 1995;
        using BoxDimType = std::array<std::array<float, 3>, 3>;

        struct FrameMeta
        {
            int32_t AtomCount;
            int32_t Step;
            float Time;
            BoxDimType BoxDim;
        };

        XtcFile() : End(Endian::current()) {}
        ~XtcFile() = default;

        XtcFile(const XtcFile&) = delete;
        XtcFile& operator=(const XtcFile&) = delete;

        void open(const char* filename);
        bool isOpen() const { return File.is_open(); }
        FrameMeta readFrameMeta();
        FrameMeta readFrame(float result[]);
        bool eof();
        void close() { File.close(); }

    private:
        std::ifstream File;
        const Endian::Endian End;

        // Return true if good, false if error or EOF
        template<typename T> bool read(T* value, size_t count=1)
        {
            if(!(File.read(reinterpret_cast<char*>(value), sizeof(T) * count)))
            {
                return false;
            }

            if(End == Endian::LITTLE)
            {
                for(size_t i = 0; i < count; i++)
                {
                    Endian::swap(*(value + i));
                }
            }
            return true;
        }

        int32_t xdrfile_decompress_coord_float(
            float* ptr, int32_t* size, float* precision);
        FrameMeta readFrameMetaAndStay();
    };

} // namespace libmd

#endif
