#ifndef MATRIX_ALGORITHMS_HPP
#define MATRIX_ALGORITHMS_HPP

#include "types.hpp"

class MatrixAlgorithms {
public:
    MatrixAlgorithms();
    Mat exp(const Mat A, int p = 20, int q = 20, double tol = 10e-12);
    Mat log(const Mat A, double tol = 10e-12);
    Vector SparseSolve(SpMat A, Vector b, double eps, int maxIter, int type = 0);
private:
    Mat padeiter_sqrt(Mat A, double tol);
    double facfrac(int i, int p, int q);
    double norm(Mat A, int p = 2);
    std::vector< std::vector<int> > pf_table;
};



#endif // MATRIX_ALGORITHMS_HPP