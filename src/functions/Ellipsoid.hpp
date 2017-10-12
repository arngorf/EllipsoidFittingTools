#ifndef ELLIPSOID_HPP
#define ELLIPSOID_HPP

#include "types.hpp"

class Ellipsoid {
public:

    Ellipsoid() {};

    Ellipsoid(Vector const &center, Vector const &radii, Mat const &evecs, Vector const &v);

    Ellipsoid(Vector const &center, Vector const &radii, Mat const &R);

    Vector getCenter() {return center;};

    double getRadii(int i) {return radii(i);};

    Vector getRadii() {return radii;};

    Mat getRotationMatrix() {return R;};

    Mat getDiagonalMatrix() {return D;};

    Mat getCoefficientMatrix() {return R.transpose()*D*R;};

    Vector getRadiiDirection(int i) {return R.row(i).transpose();};

    Vector getAlgebraicCoefficients() {return v;};

    void makeRepresentationUnique();

    void getRadii_k(double &kx1, double &ky1, double &kx2, double &ky2);

private:

    void setDiagonalMatrix(Vector const &radii);
    void sortDiagonalEntries();

    Vector center;
    Vector radii;
    Mat R;
    Mat D;
    Vector v;
};

#endif // ELLIPSOID_HPP