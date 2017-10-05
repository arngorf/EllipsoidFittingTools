#include <iostream>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "../../functions/types.hpp"
#include "../../functions/PyPlot.hpp"
#include "../../functions/MatrixAlgorithms.hpp"
#include "../../functions/MatrixHelpers.hpp"
#include "../../cuda_kernels/GpuSolver.h"

#include "Smoke2D.hpp"

using Eigen::MatrixXd;
using Eigen::MatrixBase;

/*
Matlab:             Eigen3:
*) zeros(3,2)       *) Mat::Zero(3,2)
*) y(i,j)           *) y(i-1,y-1);

*/

void printstrn(char const str[]) {
    std::cout << str << std::endl;
}

void sparseTest() {
    MatrixHelpers MatHelp;
    MatrixAlgorithms MatAlg;

    int n = 10;

    std::vector<double> values({-1,2,1});
    std::vector<int> indices({0,1,3});

    SpMat A = MatHelp.BandMatrixSparse(n, values, indices);
    Vector b(n);

    for (int i = 0; i < n; ++i) {
        b(i) = (i - n/2.0) * (i - n/2.0);
    }

    Vector x = MatAlg.SparseSolve(A, b, 10e-10, 1000);

    std::cout << A << std::endl;
    std::cout << x << std::endl;
    std::cout << A*x << std::endl;
    std::cout << b << std::endl;
}

void cudaTest() {
    GpuInterface obj;
    obj.setY(16);
    std::cout << obj.calculateSum() << std::endl;

    MatrixHelpers matHelp;
    Meshgrid G(1,3,1,3,3,3);
    Vector X = G.getAllX();
    Vector Y = G.getAllY();
    for (int i = 0; i < 9; ++i) {
        X(i) += 0;
        Y(i) += 1.0;
    }
    Mat H(3,3);
    H << 0,1,0,0.5,1,0.5,0,1,0;
    std::cout << H << std::endl;
    Mat result = matHelp.interpolate(G, H, X, Y);
    result.resize(3,3);
    std::cout << result << std::endl;
}

int main() {

    //sparseTest();
    Smoke2D sim(1, 1, 540, 540);
    sim.run();
    //cudaTest();

    return 0;
}
