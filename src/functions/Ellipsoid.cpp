#include "Ellipsoid.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <iostream>

#include "ellipsoidHelpers.hpp"


Ellipsoid::Ellipsoid(Vector const &center, Vector const &radii, Mat const &evecs, Vector const &v)
                    : center(center), radii(radii), R(evecs.transpose())
                    , D(radii.asDiagonal()), v(v)
{
    setDiagonalMatrix(radii);
    sortDiagonalEntries();
}

Ellipsoid::Ellipsoid(Vector const &center, Vector const &radii, Mat const &R)
                    : center(center), radii(radii), R(R)
                    , D(radii.asDiagonal()), v(Vector::Zero(10))
{
    setDiagonalMatrix(radii);
    sortDiagonalEntries();

    Mat H(3,3);
    H = R.transpose()*D*R;
    v(0) = H(0,0);
    v(1) = H(1,1);
    v(2) = H(2,2);
    v(3) = H(1,0);
    v(4) = H(2,0);
    v(5) = H(2,1);

    // Missing calculation of v for completion.
}

void Ellipsoid::makeRepresentationUnique()
{
    /*if (D(0,0) < D(1,1))
    {
        Mat P(3,3);
        P << 0,1,0,
             1,0,0,
             0,0,1;

        std::swap(D(0,0), D(1,1));

        R = P*R;
    }

    if (D(1,1) < D(2,2))
    {
        Mat P(3,3);
        P << 1,0,0,
             0,0,1,
             0,1,0;

        std::swap(D(1,1), D(2,2));

        R = P*R;
    }

    if (D(0,0) < D(1,1))
    {
        Mat P(3,3);
        P << 0,1,0,
             1,0,0,
             0,0,1;

        std::swap(D(0,0), D(1,1));

        R = P*R;
    }*/

    Eigen::EigenSolver<Mat> eSolve;

    Mat H = R.transpose()*D*R;

    eSolve.compute(H);
    Mat nEvecs(3,3);
    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < 3; ++i)
        {
            nEvecs(j,i) = eSolve.eigenvectors().col(j)(i).real();
        }
    }

    std::cout << "Eigenvectors:\n" << eSolve.eigenvectors() << std::endl;
    std::cout << "Eigenvalues: " << eSolve.eigenvalues()(0).real() << ", " << eSolve.eigenvalues()(1).real() << ", "  << eSolve.eigenvalues()(2).real() << ", "<< std::endl;

    std::cout << "R:\n" << R << std::endl;

    Mat Rf(3,3);
    Vector quat(4);

    double angle = PI;
    double cosRot = std::cos(angle/2);
    double sinRot = std::sin(angle/2);

    quat << cosRot, 0, sinRot, 0;// nEvecs.col(0)(0) * sinRot, nEvecs.col(0)(1) * sinRot, nEvecs.col(0)(2) * sinRot;

    Rf = QuatToMatrix(quat);
    std::cout << "Rf:\n" << Rf << std::endl;
    R = Rf*R.transpose();

    std::cout << "R:\n" << R << std::endl;
    std::cout << "H:\n" << H << std::endl;
    H = R.transpose()*D*R;
    std::cout << "H:\n" << H << std::endl;
    eSolve.compute(H);

    std::cout << "Eigenvectors:\n" << eSolve.eigenvectors() << std::endl;
    std::cout << "Eigenvalues: " << eSolve.eigenvalues()(0).real() << ", " << eSolve.eigenvalues()(1).real() << ", "  << eSolve.eigenvalues()(2).real() << ", "<< std::endl;
}

void Ellipsoid::getRadii_k(double &kx1, double &ky1, double &kx2, double &ky2)
{
    double A = v(0);
    double B = v(1);
    double C = v(2);
    double D = v(3);
    double E = v(4);
    double F = v(5);

    kx1 = - E/C;
    ky1 = - F/C;
    kx2 = (D*F - B*E) / (A*B - D*D);
    ky2 = (D*E - A*F) / (A*B - D*D);


}

void Ellipsoid::setDiagonalMatrix(Vector const &radii)
{
    D = Mat::Zero(3,3);
    for (int i = 0; i < 3; ++i)
    {
        D(i,i) = 1.0/std::pow(radii(i), 2);
    }
}

void Ellipsoid::sortDiagonalEntries()
{
    /*for (int i = 3; i > 1; --i)
    {
        for (int j = 0; j < i - 1; ++j)
        {
            if (radii(j+1) > radii(j))
            {
                R.row(j).swap(R.row(j+1));
                std::swap(radii(j),radii(j+1));
                std::swap(D(j,j),D(j+1,j+1));
            }
        }
    }*/

    /*for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (R(i,j) > 0)
            {
                break;
            }
            else if (R(i,j) < 0)
            {
                R.row(i) = -R.row(i);
                break;
            }
        }
    }*/
}