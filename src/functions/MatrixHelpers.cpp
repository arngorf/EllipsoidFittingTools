#include "MatrixHelpers.hpp"
#include <assert.h>
#include <stdexcept>
#include <iostream>

const double Meshgrid::eps = 10e-15;

Meshgrid::Meshgrid(double xBegin, double xEnd,
                   double yBegin, double yEnd,
                   int xN, int yN)
                        : xBegin(xBegin), xEnd(xEnd)
                        , yBegin(yBegin), yEnd(yEnd)
                        , xN(xN), yN(yN)
                        , dx((xEnd-xBegin)/double(xN-1))
                        , dy((yEnd-yBegin)/double(yN-1))
                        , xes(xN, yN), yes(yN, xN) {

    if (xBegin >= xEnd) throw std::invalid_argument("Meshgrid::Meshgrid:: Non-positive range in x");
    if (yBegin >= yEnd) throw std::invalid_argument("Meshgrid::Meshgrid:: Non-positive range in y");
    if (xN <= 0) throw std::invalid_argument("Meshgrid::Meshgrid:: Non-positive number of points");
    if (yN <= 0) throw std::invalid_argument("Meshgrid::Meshgrid:: Non-positive number of points");

    for (int i = 0; i < xN; ++i) {
        double newVal = (xEnd * i + xBegin * (xN - 1 - i)) / (xN - 1);
        xes(i,0) = newVal;
    }
    for (int i = 0; i < yN; ++i) {
        double newVal = (yEnd * i + yBegin * (yN - 1 - i)) / (yN - 1);
        yes(i,0) = newVal;
    }

    for (int i = 1; i < yN; ++i) {
        xes.block(0,i,xN,1) = xes.block(0,0,xN,1);
    }

    xes.transposeInPlace();

    for (int i = 1; i < xN; ++i) {
        yes.block(0,i,yN,1) = yes.block(0,0,yN,1);
    }
}

inline int Meshgrid::getPreviousX(double val) {
    int index = (val - xBegin) / dx;
    if (index < 0) return -1;
    if (index >= xN) return -1;
    return index;
}

inline int Meshgrid::getPreviousY(double val) {
    int index = (val - yBegin) / dy;
    if (index < 0) return -1;
    if (index >= yN) return -1;
    return index;
}

inline int Meshgrid::getNextX(double val) {
    return getPreviousX(val + dx);
}

inline int Meshgrid::getNextY(double val) {
    return getPreviousY(val + dy);
}

inline double Meshgrid::getX(int i) {
    if (i < 0 || i >= xN) {
        throw std::out_of_range("Meshgrid::getX: Index out of range");
    }
    return xes(0,i);
}

inline double Meshgrid::getY(int i) {
    if (i < 0 || i >= yN) {
        throw std::out_of_range("Meshgrid::getY: Index out of range");
    }
    return yes(i,0);
}

Vector Meshgrid::getX() {
    Mat result = xes.block(0, 0, 1, xN);
    result.transposeInPlace();
    return result;
}

Vector Meshgrid::getY() {
    Vector result = yes.block(0, 0, yN, 1);
    return result;
}

Mat Meshgrid::getAllX() {
    return xes;
}

Mat Meshgrid::getAllY() {
    return yes;
}

SpMat MatrixHelpers::BandMatrixSparse(int n, std::vector<double> values, std::vector<int> indices) {

    std::vector<Triplet> coefficients;
    for (int indexIndex = 0; indexIndex < indices.size(); ++indexIndex) {

        int index = indices[indexIndex];
        for (int i = 0; i < n - index; ++i) {
            coefficients.push_back(Triplet(index + i, i, values[indexIndex]));
            if (index != 0) {
                coefficients.push_back(Triplet(i, index + i, values[indexIndex]));
            }
        }
    }

    SpMat A(n,n);
    A.setFromTriplets(coefficients.begin(), coefficients.end());
    return A;
}

Mat MatrixHelpers::interpolate(Meshgrid G, Mat H, Vector nx, Vector ny, double fill) {

    Mat result(nx.size(), 1); // output result

    // Prev values left/right of top/bottom of wanted point.

    double dx = G.dx;
    double dy = G.dy;

    for (int i = 0; i < nx.size(); ++i) {

        int xPrevIndex = G.getPreviousX(nx[i]);
        int xNextIndex = G.getNextX(nx[i]);
        int yPrevIndex = G.getPreviousY(ny[i]);
        int yNextIndex = G.getNextY(ny[i]);

        // Completely outside grid
        if ((xPrevIndex == -1 && xNextIndex == -1) ||
            (yPrevIndex == -1 && yNextIndex == -1)) {
            result(i, 0) = fill;
            continue;
        }

        double xAlpha, yAlpha;
        double leX, leY;
        double TLH, TRH, BLH, BRH;

        // Edge cases
        if (xPrevIndex == -1) {

            // Case: Just outside left x edge
            leX = G.getX(0) - dx;
            TLH = fill;
            BLH = fill;

            if (yPrevIndex == -1) {

                // Case: Just outside top y edge - Corner case
                leY = G.getY(0) - dy;
                TRH = fill;
                BRH = H(yNextIndex, xNextIndex);

            } else if (yNextIndex == -1) {

                // Case: Just outside bottom y edge - Corner case
                leY = G.getY(yPrevIndex);
                TRH = H(yPrevIndex, xNextIndex);
                BRH = fill;

            } else {

                // Case: Not at top or bottom edge
                leY = G.getY(yPrevIndex);
                BRH = H(yNextIndex, xNextIndex);
                TRH = H(yPrevIndex, xNextIndex);

            }
        } else if (xNextIndex == -1) {

            // Case: Just outside right x edge
            leX = G.getX(xPrevIndex);
            TRH = fill;
            BRH = fill;

            if (yPrevIndex == -1) {

                // Case: Just outside top y edge - Corner case
                leY = G.getY(0) - dy;
                TLH = fill;
                BLH = H(yNextIndex, xPrevIndex);

            } else if (yNextIndex == -1) {

                // Case: Just outside bottom y edge - Corner case
                leY = G.getY(yPrevIndex);
                TLH = H(yPrevIndex, xPrevIndex);
                BLH = fill;

            } else {

                // Case: Not at top or bottom edge
                leY = G.getY(yPrevIndex);
                BLH = H(yNextIndex, xPrevIndex);
                TLH = H(yPrevIndex, xPrevIndex);

            }
        } else {

            // Case: Not at right or left edge
            leX = G.getX(xPrevIndex);

            if (yPrevIndex == -1) {

                // Case: Just outside top y edge
                leY = G.getY(0) - dy;
                TLH = fill;
                TRH = fill;
                BLH = H(yNextIndex, xPrevIndex);
                BRH = H(yNextIndex, xNextIndex);

            } else if (yNextIndex == -1) {

                // Case: Just outside bottom y edge
                leY = G.getY(yPrevIndex);
                TLH = H(yPrevIndex, xPrevIndex);
                TRH = H(yPrevIndex, xNextIndex);
                BLH = fill;
                BRH = fill;

            } else {

                // Case: Center case
                leY = G.getY(yPrevIndex);
                TLH = H(yPrevIndex, xPrevIndex);
                TRH = H(yPrevIndex, xNextIndex);
                BLH = H(yNextIndex, xPrevIndex);
                BRH = H(yNextIndex, xNextIndex);

            }

        }

        xAlpha = (nx[i] - leX) / dx;
        yAlpha = (ny[i] - leY) / dy;

        result(i, 0) = (1 - yAlpha) * (1 - xAlpha) * TLH
                     + (1 - yAlpha) *      xAlpha  * TRH
                     +      yAlpha  * (1 - xAlpha) * BLH
                     +      yAlpha  *      xAlpha  * BRH;
    }
    return result;
}