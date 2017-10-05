#include "TiffReader.hpp"

#include <iostream>
#include <sstream>

void DummyHandler(const char* module, const char* fmt, va_list ap) {
    // ignore errors and warnings (or handle them your own way)
}

int TiffReader::read(char const path[], double *image,
                     int _xOffset, int _yOffset, int _zOffset,
                     int _width, int _height, int _depth) {

    // disable warnings
    TIFFSetWarningHandler(DummyHandler);


    TIFF *tif = TIFFOpen(path, "r");
    if (tif) {

        unsigned int width, height;
        // get the size of the tiff

        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);

        if (_xOffset + _width > width) {
            std::cout << "width of image exceeded" << std::endl;
            return -1;
        }

        if (_yOffset + _height > height) {
            std::cout << "height of image exceeded" << std::endl;
            return -1;
        }

        size_t pageCounter = 0;

        do {

            uint32* raster;
            pageCounter++;
            //std::cout << pageCounter << std::endl;

            uint npixels = width*height; // get the total number of pixels

            raster = (uint32*)_TIFFmalloc(npixels * sizeof(uint32)); // allocate temp memory (must use the tiff library malloc)

            if (raster == NULL) {
                TIFFClose(tif);
                std::cout << "Could not allocate memory for raster of TIFF image\n";
                return -1;
            }

            // Check the tif read to the raster correctly
            if (!TIFFReadRGBAImage(tif, width, height, raster, 0)) {

                TIFFClose(tif);
                std::cout << "Could not read raster of TIFF image\n";
                 _TIFFfree(raster);
                free(image);
                return -1;
            }

            uint subImageSize = _height * _width;

            // iterate through all the pixels of the tif

            for (uint x = 0; x < width; x++) {
                for (uint y = 0; y < height; y++) {

                    int subImageZ = pageCounter - _zOffset;
                    int subImageY = (height - y) - _yOffset;
                    int subImageX = x - _xOffset;

                    if (subImageZ >= 0 and subImageZ < _depth  and
                        subImageY >= 0 and subImageY < _height and
                        subImageX >= 0 and subImageX < _width) {

                        // read the current pixel of the TIF
                        uint32& TiffPixel = raster[y*width + x];

                        // Take last 8 bits
                        unsigned low8bits = TiffPixel & 0xFF;

                        image[subImageSize*subImageZ + subImageY*_width + subImageX] = (double) low8bits;
                    }
                }
            }
            _TIFFfree(raster); // release temp memory

        } while (TIFFReadDirectory(tif)); // get the next tif

    if (_zOffset + _depth > pageCounter) {
        std::cout << "depth of image exceeded" << std::endl;
        free(image);
        return -1;
    }

    TIFFClose(tif);

    _depth = pageCounter;
    _height = height;
    _width = width;


    } else {
        std::cout << "Error reading file " << path << std::endl;
        return -1;
    }

    return 0;

}