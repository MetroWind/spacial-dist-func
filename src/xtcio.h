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
