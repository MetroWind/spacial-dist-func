#include <exception>
#include <vector>

#include "xtcio.h"

namespace libmd
{
    namespace                   // Ugly stuff
    {

        // Shamelessly copied from
        // https://github.com/wesbarnett/libxdrfile/blob/master/src/xdrfile.c.
        //
        /* Internal support routines for reading/writing compressed
         * coordinates sizeofint - calculate smallest number of bits
         * necessary to represent a certain integer.
         */
        static uint32_t sizeofint(int32_t size)
        {
            uint32_t num = 1;
            uint32_t num_of_bits = 0;

            while(static_cast<uint32_t>(size) >= num && num_of_bits < 32)
            {
                num_of_bits++;
                num <<= 1;
            }
            return num_of_bits;
        }


        // Shamelessly copied from
        // https://github.com/wesbarnett/libxdrfile/blob/master/src/xdrfile.c.
        //
/*
 * sizeofints - calculate 'bitsize' of compressed ints
 *
 * given a number of small unsigned integers and the maximum value
 * return the number of bits needed to read or write them with the
 * routines encodeints/decodeints. You need this parameter when
 * calling those routines.
 * (However, in some cases we can just use the variable 'smallidx'
 * which is the exact number of bits, and them we dont need to call
 * this routine).
 */
        static uint32_t sizeofints(uint32_t num_of_ints, uint32_t sizes[])
        {
            uint32_t i, num;
            uint32_t num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
            num_of_bytes = 1;
            bytes[0] = 1;
            num_of_bits = 0;
            for (i=0; i < num_of_ints; i++)
            {
                tmp = 0;
                for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++)
                {
                    tmp = bytes[bytecnt] * sizes[i] + tmp;
                    bytes[bytecnt] = tmp & 0xff;
                    tmp >>= 8;
                }
                while (tmp != 0)
                {
                    bytes[bytecnt++] = tmp & 0xff;
                    tmp >>= 8;
                }
                num_of_bytes = bytecnt;
            }
            num = 1;
            num_of_bytes--;
            while (bytes[num_of_bytes] >= num)
            {
                num_of_bits++;
                num *= 2;
            }
            return num_of_bits + num_of_bytes * 8;

        }

        // Shamelessly copied from
        // https://github.com/wesbarnett/libxdrfile/blob/master/src/xdrfile.c.
        //
/*
 * decodebits - decode number from buf using specified number of bits
 *
 * extract the number of bits from the array buf and construct an integer
 * from it. Return that value.
 *
 */

        static int32_t decodebits(int32_t buf[], int32_t num_of_bits)
        {
            int32_t cnt, num;
            uint32_t lastbits, lastbyte;
            unsigned char * cbuf;
            int32_t mask = (1 << num_of_bits) -1;
            cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
            cnt = buf[0];
            lastbits = (uint32_t) buf[1];
            lastbyte = (uint32_t) buf[2];

            num = 0;
            while (num_of_bits >= 8)
            {
                lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
                num |=  (lastbyte >> lastbits) << (num_of_bits - 8);
                num_of_bits -=8;
            }
            if (num_of_bits > 0)
            {
                if (lastbits < static_cast<uint32_t>(num_of_bits))
                {
                    lastbits += 8;
                    lastbyte = (lastbyte << 8) | cbuf[cnt++];
                }
                lastbits -= num_of_bits;
                num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
            }
            num &= mask;
            buf[0] = cnt;
            buf[1] = lastbits;
            buf[2] = lastbyte;
            return num;
        }

        // Shamelessly copied from
        // https://github.com/wesbarnett/libxdrfile/blob/master/src/xdrfile.c.
        //
/*
 * decodeints - decode 'small' integers from the buf array
 *
 * this routine is the inverse from encodeints() and decodes the small integers
 * written to buf by calculating the remainder and doing divisions with
 * the given sizes[]. You need to specify the total number of bits to be
 * used from buf in num_of_bits.
 *
 */

        static void decodeints(int32_t buf[], int32_t num_of_ints, int32_t num_of_bits,
                   uint32_t sizes[], int32_t nums[])
        {
            int32_t bytes[32];
            int32_t i, j, num_of_bytes, p, num;

            bytes[1] = bytes[2] = bytes[3] = 0;
            num_of_bytes = 0;
            while (num_of_bits > 8)
            {
                int32_t x = decodebits(buf, 8);
                bytes[num_of_bytes++] = x;
                num_of_bits -= 8;
            }
            if (num_of_bits > 0)
            {
                bytes[num_of_bytes++] = decodebits(buf, num_of_bits);
            }
            for (i = num_of_ints-1; i > 0; i--)
            {
                num = 0;
                for (j = num_of_bytes-1; j >=0; j--)
                {
                    num = (num << 8) | bytes[j];
                    p = num / sizes[i];
                    bytes[j] = p;
                    num = num - p * sizes[i];
                }
                nums[i] = num;
            }
            nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
        }
    } // namespace

    // Shamelessly copied from
    // https://github.com/wesbarnett/libxdrfile/blob/master/src/xdrfile.c.
    //
    /* Compressed coordinate routines - modified from the original
     * implementation by Frans v. Hoesel to make them threadsafe.
     */
    int32_t XtcFile :: xdrfile_decompress_coord_float(
        float* ptr, int* size, float* precision)
    {
        if(ptr == nullptr)
            return -1;

        constexpr int32_t magicints[] =
            {
                0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
                80, 101, 128, 161, 203, 256, 322, 406, 512, 645, 812, 1024, 1290,
                1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321, 13003,
                16384, 20642, 26007, 32768, 41285, 52015, 65536,82570, 104031,
                131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
                832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021,
                4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216
            };

        constexpr int32_t FIRSTIDX = 9;

        /* note that magicints[FIRSTIDX-1] == 0 */
        // constexpr int32_t LASTIDX = (sizeof(magicints) / sizeof(*magicints));

        int32_t minint[3], maxint[3], *lip;
        int32_t smallidx;
        uint32_t sizeint[3], sizesmall[3], bitsizeint[3], size3;
        int32_t k, *buf1, lsize, flag;
        int32_t smallnum, smaller, i, is_smaller, run;
        float *lfp, inv_precision;
        int32_t tmp, *thiscoord,  prevcoord[3];
        uint32_t bitsize;

        bitsizeint[0] = 0;
        bitsizeint[1] = 0;
        bitsizeint[2] = 0;

        read(&lsize);
        if (*size < lsize)
        {
            fprintf(stderr, "Requested to decompress %d coords, file contains %d\n",
                    *size, lsize);
            return -1;
        }
        *size = lsize;
        size3 = *size * 3;

        std::vector<int> Buf1Body(size3, 0);
        buf1 = Buf1Body.data();
        // buf2 = (int *)malloc(sizeof(int)*buf2size);
        std::array<int, 3> Buf2Header;

            /* Dont bother with compression for three atoms or less */
            if(*size<=9)
            {
                read(ptr, size3);
                /* return number of coords, not floats */
                return size3 / 3;
            }

        /* Compression-time if we got here. Read precision first */
        read(precision);

        /* Buf2Header[0-2] are special and do not contain actual data */
        Buf2Header[0] = Buf2Header[1] = Buf2Header[2] = 0;
        read(minint, 3);
        read(maxint, 3);

        sizeint[0] = maxint[0] - minint[0]+1;
        sizeint[1] = maxint[1] - minint[1]+1;
        sizeint[2] = maxint[2] - minint[2]+1;

        /* check if one of the sizes is to big to be multiplied */
        if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff)
        {
            bitsizeint[0] = sizeofint(sizeint[0]);
            bitsizeint[1] = sizeofint(sizeint[1]);
            bitsizeint[2] = sizeofint(sizeint[2]);
            bitsize = 0; /* flag the use of large sizes */
        }
        else
        {
            bitsize = sizeofints(3, sizeint);
        }

        read(&smallidx);

        tmp = smallidx+8;
        tmp = smallidx-1;
        tmp = (FIRSTIDX>tmp) ? FIRSTIDX : tmp;
        smaller = magicints[tmp] / 2;
        smallnum = magicints[smallidx] / 2;
        sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;

        /* Buf2Header[0] holds the length in bytes */

        read(Buf2Header.data());
        // Not in original xdrfile.c: Pad to 4 bytes...?
        if(Buf2Header[0] % 4 != 0)
        {
            Buf2Header[0] = (Buf2Header[0] / 4 + 1) * 4;
        }

        std::vector<int> Buf2Body(Buf2Header[0]/sizeof(int32_t) + 3 + 1, 0);
        int* buf2 = Buf2Body.data();

        File.read(reinterpret_cast<char*>(&(buf2[3])),
                  static_cast<uint32_t>(Buf2Header[0]));
        // if (xdrfile_read_opaque((char *)&(buf2[3]),(unsigned int)buf2[0],xfp) == 0)
        //     return 0;
        // buf2[0] = buf2[1] = buf2[2] = 0;

        lfp = ptr;
        inv_precision = 1.0 / * precision;
        run = 0;
        i = 0;
        lip = buf1;
        while ( i < lsize )
        {
            thiscoord = (int32_t *)(lip) + i * 3;

            if (bitsize == 0)
            {
                thiscoord[0] = decodebits(buf2, bitsizeint[0]);
                thiscoord[1] = decodebits(buf2, bitsizeint[1]);
                thiscoord[2] = decodebits(buf2, bitsizeint[2]);
            }
            else
            {
                decodeints(buf2, 3, bitsize, sizeint, thiscoord);
            }
            i++;
            thiscoord[0] += minint[0];
            thiscoord[1] += minint[1];
            thiscoord[2] += minint[2];

            prevcoord[0] = thiscoord[0];
            prevcoord[1] = thiscoord[1];
            prevcoord[2] = thiscoord[2];

            flag = decodebits(buf2, 1);
            is_smaller = 0;
            if (flag == 1)
            {
                run = decodebits(buf2, 5);
                is_smaller = run % 3;
                run -= is_smaller;
                is_smaller--;
            }
            if (run > 0)
            {
                thiscoord += 3;
                for (k = 0; k < run; k+=3)
                {
                    decodeints(buf2, 3, smallidx, sizesmall, thiscoord);
                    i++;
                    thiscoord[0] += prevcoord[0] - smallnum;
                    thiscoord[1] += prevcoord[1] - smallnum;
                    thiscoord[2] += prevcoord[2] - smallnum;
                    if (k == 0) {
                        /* interchange first with second atom for better
                         * compression of water molecules
                         */
                        tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
                        prevcoord[0] = tmp;
                        tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
                        prevcoord[1] = tmp;
                        tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
                        prevcoord[2] = tmp;
                        *lfp++ = prevcoord[0] * inv_precision;
                        *lfp++ = prevcoord[1] * inv_precision;
                        *lfp++ = prevcoord[2] * inv_precision;
                    } else {
                        prevcoord[0] = thiscoord[0];
                        prevcoord[1] = thiscoord[1];
                        prevcoord[2] = thiscoord[2];
                    }
                    *lfp++ = thiscoord[0] * inv_precision;
                    *lfp++ = thiscoord[1] * inv_precision;
                    *lfp++ = thiscoord[2] * inv_precision;
                }
            }
            else
            {
                *lfp++ = thiscoord[0] * inv_precision;
                *lfp++ = thiscoord[1] * inv_precision;
                *lfp++ = thiscoord[2] * inv_precision;
            }
            smallidx += is_smaller;
            if (is_smaller < 0)
            {
                smallnum = smaller;

                if (smallidx > FIRSTIDX)
                {
                    smaller = magicints[smallidx - 1] /2;
                }
                else
                {
                    smaller = 0;
                }
            }
            else if (is_smaller > 0)
            {
                smaller = smallnum;
                smallnum = magicints[smallidx] / 2;
            }
            sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
        }
        return *size;
    }

    void XtcFile :: open(const char* filename)
    {
        File.open(filename);
    }

    bool XtcFile :: eof()
    {
        return File.peek() == std::char_traits<char>::eof();
    }

    XtcFile::FrameMeta XtcFile :: readFrameMetaAndStay()
    {
        int32_t magic;
        read(&magic);
        if(magic != MAGIC)
        {
            throw std::runtime_error("incorrect magic");
        }

        XtcFile::FrameMeta Meta;
        read(&Meta.AtomCount);
        read(&Meta.Step);
        read(&Meta.Time);
        for(size_t i = 0; i < 3; i++)
        {
            for(size_t j = 0; j < 3; j++)
            {
                read(&(Meta.BoxDim[i][j]));
            }
        }
        return Meta;
    }

    XtcFile::FrameMeta XtcFile :: readFrameMeta()
    {
        auto FrameBegin = File.tellg();
        auto Meta = readFrameMetaAndStay();
        File.seekg(FrameBegin);
        return Meta;
    }

    XtcFile::FrameMeta XtcFile :: readFrame(float result[])
    {
        auto Meta = readFrameMetaAndStay();
        float precision;
        xdrfile_decompress_coord_float(result, &Meta.AtomCount, &precision);

        return Meta;
    }

}
