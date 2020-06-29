#ifndef SDF_ENDIAN_H
#define SDF_ENDIAN_H

#include <algorithm>

namespace libmd
{
    namespace Endian
    {
        enum Endian { BIG, LITTLE };
        inline Endian current()
        {
            union{ uint16_t i; char c[2]; } bint = {0x0102};

            if(bint.c[0] == 1)
            {
                return BIG;
            }
            else
            {
                return LITTLE;
            }
        }

        template <typename T> T swapped(const T& x)
        {
            union U
            {
                T val;
                std::array<std::uint8_t, sizeof(T)> raw;
            } src, dst;

            src.val = x;
            std::reverse_copy(src.raw.begin(), src.raw.end(), dst.raw.begin());
            return dst.val;
        }

        template <typename T> void swap(T& x)
        {
            union U
            {
                T val;
                std::array<std::uint8_t, sizeof(T)> raw;
            } src, dst;
            src.val = x;
            std::reverse_copy(src.raw.begin(), src.raw.end(), dst.raw.begin());
            x = dst.val;
        }
    }
}

#endif
