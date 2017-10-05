#ifndef TIFFREADER_HPP
#define TIFFREADER_HPP

#include <tiffio.h>
#include "types.hpp"

void DummyHandler(const char* module, const char* fmt, va_list ap);

class TiffReader {

public:

    int read(char const path[], double *image,
             int _xOffset, int _yOffset, int _zOffset,
             int _width, int _height, int _depth);
};

#endif // TIFFREADER_HPP