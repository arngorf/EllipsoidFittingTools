#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

//#include "DistPointHyperellipsoid.hpp"
#include "ellipsoidHelpers.hpp"
#include "EllipsoidLeastSquareFit.hpp"
#include "EllipsoidMinimizer.hpp"
//#include "MatrixAlgorithms.hpp"
//#include "MatrixHelpers.hpp"
#include "Minimizer.hpp"
#include "MinimizerN.hpp"
#include "PyPlot.hpp"

#include "ChanVese3D.hpp"
//#include "Rotation.hpp"
//#include "TaskMaster.hpp"
#include "TiffReader.hpp"
#include "types.hpp"
#include <stdlib.h>

//using Eigen::MatrixXd;
//using Eigen::MatrixBase;

void testChanVese() {

    TiffReader tr;
    ChanVese3D cv3d;
    Graph2D *plt = Graph2D::Instance();

    //size_t xVesicleCenter = 1197, yVesicleCenter = 863, zVesicleCenter = 165;
    size_t xVesicleCenter = 1038, yVesicleCenter = 555, zVesicleCenter = 362;
    size_t width = 24, height = 20, depth = 24;
    size_t xOffset = xVesicleCenter - width/2, yOffset = yVesicleCenter - height/2, zOffset = zVesicleCenter - depth/2;
    size_t arrayLength = width * height * depth;
    size_t paddedArrayLength = (width + 2) * (height + 2) * (depth + 2);
    double *image = (double*) malloc(arrayLength * sizeof(double));
    double *phaseField = (double*) malloc(paddedArrayLength * sizeof(double));

    int readErr = tr.read("/home/dith/MahdiehVesicle/3dvesicleSnakes/LaussaneDataNAnnotations/volumedata.tif",
                          image, xOffset, yOffset, zOffset, width, height, depth);

    if (readErr != 0) {
        std::cout << "Error reading image\n";
        exit(0);
    }

    Mat imSlice(height,width);
    Mat phiSlice(height,width);
    std::cout << height << ", " << width << std::endl;

    cv3d.segmentImageCPU(image, phaseField, depth, height, width);

    double ppMin = 0;
    double ppMax = 0;

    for (int z = 0; z < depth; ++z) {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                imSlice(y, x) = image[gIndex(z,y,x,height,width,0)];
                phiSlice(y, x) = phaseField[gIndex(z,y,x,height,width,1)];
                ppMin = std::min(ppMin, phaseField[gIndex(z,y,x,height,width,1)]);
                ppMax = std::max(ppMax, phaseField[gIndex(z,y,x,height,width,1)]);
            }
        }

        plt->subplot(4, 6, z + 1);
        plt->title("slice " + STR(zOffset + z));
        plt->imshow(imSlice, phiSlice, 0, 255);
    }
    plt->show();

    for (int z = 0; z < depth; ++z) {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                phiSlice(y, x) = phaseField[gIndex(z,y,x,height,width,1)];

            }
        }

        plt->subplot(4, 6, z + 1);
        plt->title("slice " + STR(zOffset + z));
        plt->imshow(phiSlice, ppMin, ppMax);
    }
    plt->show();

    Mat X(2*(width*height + width*depth + height*depth), 3);

    int curIndex = 0;
    for (int z = 1; z < depth; ++z) {
        for (int y = 1; y < height; ++y) {
            for (int x = 1; x < width; ++x) {

                double b = phaseField[gIndex(z,y,x,height,width,1)];
                double a = phaseField[gIndex(z,y,x-1,height,width,1)];

                if (a * b < 0) {
                    double interpolant = a/(a-b);
                    X(curIndex, 0) = x - 1 + interpolant;
                    X(curIndex, 1) = y;
                    X(curIndex, 2) = z;
                    ++curIndex;
                }

                a = phaseField[gIndex(z,y-1,x,height,width,1)];

                if (a * b < 0) {
                    double interpolant = a/(a-b);
                    X(curIndex, 0) = x;
                    X(curIndex, 1) = y - 1 + interpolant;
                    X(curIndex, 2) = z;
                    ++curIndex;
                }

                a = phaseField[gIndex(z-1,y,x,height,width,1)];

                if (a * b < 0) {
                    double interpolant = a/(a-b);
                    X(curIndex, 0) = x;
                    X(curIndex, 1) = y;
                    X(curIndex, 2) = z - 1 + interpolant;
                    ++curIndex;
                }
            }
        }
    }

    plt->scatter3d(X.block(0,0,curIndex,3));
    plt->show();

    free(image);
    free(phaseField);
}

int main() {

    testChanVese();

    return 0;
}