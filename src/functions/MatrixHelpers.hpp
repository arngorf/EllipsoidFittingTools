#ifndef MATRIX_HELPERS_HPP
#define MATRIX_HELPERS_HPP

#include "types.hpp"

class Meshgrid {
public:

    Meshgrid(double xBegin, double xEnd,
             double yBegin, double yEnd,
             int xN, int yN);

    int getPreviousX(double val);
    int getPreviousY(double val);
    int getNextX(double val);
    int getNextY(double val);

    double getX(int i);
    double getY(int i);

    Vector getX();
    Vector getY();
    Mat getAllX();
    Mat getAllY();

    const double xBegin;
    const double xEnd;
    const double yBegin;
    const double yEnd;
    const int xN;
    const int yN;
    const double dx;
    const double dy;

private:

    Mat xes;
    Mat yes;

    static const double eps;
};

class MatrixHelpers {
public:

    MatrixHelpers() {};

    SpMat BandMatrixSparse(int n, std::vector<double> values, std::vector<int> indices);

    Mat interpolate(Meshgrid G, Mat H, Vector nx, Vector ny, double fill = 0.0);

};

#endif // MATRIX_HELPERS_HPP