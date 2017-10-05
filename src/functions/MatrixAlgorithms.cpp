
#include "MatrixAlgorithms.hpp"
#include <iostream>
#include <algorithm>
#include <exception>
#include <Eigen/Dense>

MatrixAlgorithms::MatrixAlgorithms() {
    // Prime factors
    // You do crazy things for a little bit of speed :o)
    pf_table.push_back({0});           pf_table.push_back({});      // 0, 1
    pf_table.push_back({2});           pf_table.push_back({3});     // 2, 3
    pf_table.push_back({2,2});         pf_table.push_back({5});     // 4, 5
    pf_table.push_back({2,3});         pf_table.push_back({7});     // 6, 7
    pf_table.push_back({2,2,2});       pf_table.push_back({3,3});   // 8, 9
    pf_table.push_back({2,5});         pf_table.push_back({11});    // 10, 11
    pf_table.push_back({2,2,3});       pf_table.push_back({13});    // 12, 13
    pf_table.push_back({2,7});         pf_table.push_back({3,5});   // 14, 15
    pf_table.push_back({2,2,2,2});     pf_table.push_back({17});    // 16, 17
    pf_table.push_back({2,3,3});       pf_table.push_back({19});    // 18, 19
    pf_table.push_back({2,2,5});       pf_table.push_back({3,7});   // 20, 21
    pf_table.push_back({2,11});        pf_table.push_back({23});    // 22, 23
    pf_table.push_back({2,2,2,3});     pf_table.push_back({5,5});   // 24, 25
    pf_table.push_back({2,13});        pf_table.push_back({3,3,3}); // 26, 27
    pf_table.push_back({2,2,7});       pf_table.push_back({29});    // 28, 29
    pf_table.push_back({2,3,5});       pf_table.push_back({31});    // 30, 31
    pf_table.push_back({2,2,2,2,2});   pf_table.push_back({3,11});  // 32, 33
    pf_table.push_back({2,17});        pf_table.push_back({5,7});   // 34, 35
    pf_table.push_back({2,2,3,3});     pf_table.push_back({37});    // 36, 37
    pf_table.push_back({2,19});        pf_table.push_back({3,13});  // 38, 39
    pf_table.push_back({2,2,2,5});                                  // 40
}

Mat MatrixAlgorithms::exp(const Mat A_in, int p, int q, double tol) {

    if (A_in.cols() != A_in.rows()) {
        throw std::invalid_argument("Matrix must be square");
    }

    int n = A_in.rows();

    double d = norm(A_in, 1);

    int k = 0;

    while (d > tol) {
        d /= 2;
        ++k;
    }

    d = 1.0/double(pow(2,k));

    // Scaling
    Mat A = A_in*d;

    Mat N = Mat::Identity(n,n);
    Mat D = Mat::Identity(n,n);

    Mat An = Mat::Identity(n,n);
    for (int i = 1; i < p; ++i) {
        An *= A;
        N += facfrac(i,p,q) * An;
    }

    An = Mat::Identity(n,n);
    for (int j = 1; j < q; ++j) {
        An *= -A;
        D += facfrac(j,q,p) * An;
    }

    Eigen::FullPivLU<Mat> LU(D);

    Mat Exp = LU.solve(N);

    // Squaring
    for (int i = 0; i < k; ++i) {
        Exp *= Exp;
    }

    return Exp;

    //Mat L = LU.matrixLU();
    //Mat U = LU.matrixLU().triangularView<Upper>();
    //Mat P = Mat::Identity(n,n);
    //Mat Q = Mat::Identity(n,n);

    //L.block<5,3>(0,0).triangularView<StrictlyLower>() = lu.matrixLU();
    //cout << l << endl;
    //cout << "Here is the U part:" << endl;
    //Matrix5x3 u = ;
    //cout << u << endl;
    //cout << "Let us now reconstruct the original matrix m:" << endl;
    //cout << lu.permutationP().inverse() * l * u * lu.permutationQ().inverse() << endl;

}

Mat MatrixAlgorithms::log(const Mat A, double tol) {
    /*Mat X = newton_sqrt(A, tol);
    if (norm(X*X - A, 2) < tol) {
        return X;
    }*/

    int k = 15;
    int n = A.rows();

    Mat X = A;

    for (int i = 0; i < k; ++i) {
        X = padeiter_sqrt(X, tol);
    }

    X -= Mat::Identity(n,n);

    X *= pow(2,k);

    //std::cout << "MatrixAlgorithms::log: Warning, no method gave results wit required tolerance level." << std::endl;

    //X = -X;// + Mat::Identity(n,n);

    return X;
}

Vector MatrixAlgorithms::SparseSolve(SpMat A, Vector b, double eps, int maxIter, int type) {
    Vector x;

    if (type == 0) {
        Eigen::SimplicialCholesky<SpMat> chol(A);
        Vector x = chol.solve(b);
    } else if (type == 1) {
        Eigen::BiCGSTAB<SpMat> solver;
        solver.setTolerance(eps);
        solver.setMaxIterations(maxIter);
        solver.compute(A);
        x = solver.solve(b);
    } else if (type == 2) {
        /*Eigen::LeastSquaresConjugateGradient<SpMat> lscg;
        lscg.setTolerance(eps);
        lscg.setMaxIterations(maxIter);
        lscg.compute(A);
        x = lscg.solve(b);*/
    }

    return x;
}

Mat MatrixAlgorithms::padeiter_sqrt(Mat A, double tol) {
    int p = 10;
    int n = A.rows();
    Mat Y = A;
    Mat Z = Mat::Identity(n,n);
    Mat oldY;
    Mat oldZ;

    std::vector<double> xi;
    std::vector<double> a;

    for (int i = 1; i <= p; ++i) {
        xi.push_back((1 + std::cos((2*i - 1)*PI/(2.0*p)))/2.0);
        a.push_back(1.0/xi[i-1] - 1);
    }

    for (int i = 0; i < 500; ++i) {
        if (/*i % 10 == 0 &&*/ norm(Y*Y-A, 2) < tol) {
            std::cout << "MatrixAlgorithms::padeiter_sqrt: Stopped at iteration " << i << " Error: " << norm(Y*Y-A, 2) << std::endl;
            return Y;
        }

        oldY = Y;
        oldZ = Z;
        Y.setZero(n,n);
        Z.setZero(n,n);

        double detY = oldY.determinant();
        double detZ = oldZ.determinant();

        if (detY == 0) {
            std::cout << "Halting: Y determinant is zero" << std::endl;
            break;
        }

        if (detZ == 0) {
            std::cout << "Halting: Z determinant is zero" << std::endl;
            break;
        }

        double g = std::abs(pow(detY * detZ, -1.0/(2.0*n)));

        for (int i = 0; i < p; ++i) {
            Y += 1.0/xi[i] * (pow(g,2)*oldZ*oldY + a[i]*Mat::Identity(n,n)).inverse();
            Z += 1.0/xi[i] * (pow(g,2)*oldY*oldZ + a[i]*Mat::Identity(n,n)).inverse();
        }

        Y = g / p * oldY * Y;
        Z = g / p * oldZ * Z;

    }

    std::cout << "MatrixAlgorithms::padeiter_sqrt: Warning (2 norm) error at: " << norm(Y*Y-A, 2) << std::endl;

    return Y;

}

double MatrixAlgorithms::facfrac(int i, int p, int q) {

    if (p + q > 40) {
        throw std::invalid_argument("p + q must be at most 40.");
    }

    if (p < i) {
        throw std::invalid_argument("i must be at most p.");
    }

    std::vector<int> numer_facts({1});
    std::vector<int> denom_facts({1});

    // numerator  : (p+q-i)!
    // denominator: (p+q)!
    for (int n = std::max(2,p+q-i+1); n <= p+q; ++n) {
        denom_facts.insert(denom_facts.end(), pf_table[n].begin(), pf_table[n].end());
    }

    // numerator  : p!
    // denominator: (p-i)!
    for (int n = std::max(p-i+1,2); n <= p; ++n) {
        numer_facts.insert(numer_facts.end(), pf_table[n].begin(), pf_table[n].end());
    }

    // denominator: i!
    for (int n = 2; n <= i; ++n) {
        denom_facts.insert(denom_facts.end(), pf_table[n].begin(), pf_table[n].end());
    }

    // Cancellation
    std::sort(numer_facts.begin(), numer_facts.end());
    std::sort(denom_facts.begin(), denom_facts.end());

    int n = 0;
    int d = 0;
    while (n < numer_facts.size() && d < denom_facts.size()) {

        if (numer_facts[n] == 1) {
            ++n;
            continue;
        }
        if (denom_facts[d] == 1) {
            ++d;
            continue;
        }
        if (numer_facts[n] == denom_facts[d]) {
            numer_facts[n] = 1;
            denom_facts[d] = 1;
        }
        if (numer_facts[n] < denom_facts[d]) {
            ++n;
        } else {
            ++d;
        }
    }

    // Multiply and divide
    double result = 1.0;
    n = 0;
    d = 0;
    while (n < numer_facts.size() || d < denom_facts.size()) {
        if ((result < 100 && n < numer_facts.size()) || d == denom_facts.size()) {
            result *= numer_facts[n];
            ++n;
        } else {
            result /= denom_facts[d];
            ++d;
        }
    }

    return result;
}

double MatrixAlgorithms::norm(Mat A, int p) {

    double d = 0;

    for (int j = 0; j < A.rows(); ++j) {
        for (int i = 0; i < A.cols(); ++i)  {
            d += pow(std::abs(A(j,i)), p);
        }
    }

    return pow(d,1.0/p);
}
