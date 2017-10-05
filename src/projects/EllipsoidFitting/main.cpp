#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <mutex>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "../../cuda_kernels/GpuSolver.h"
#include "../../functions/DistPointHyperellipsoid.hpp"
#include "../../functions/ellipsoidHelpers.hpp"
#include "../../functions/EllipsoidLeastSquareFit.hpp"
#include "../../functions/EllipsoidMinimizer.hpp"
#include "../../functions/MatrixAlgorithms.hpp"
#include "../../functions/MatrixHelpers.hpp"
#include "../../functions/Minimizer.hpp"
#include "../../functions/MinimizerN.hpp"
#include "../../functions/PyPlot.hpp"
#include "../../functions/Rotation.hpp"
#include "../../functions/types.hpp"
#include "../../functions/TaskMaster.hpp"

using Eigen::MatrixXd;
using Eigen::MatrixBase;

void printstrn(char const str[]) {
    std::cout << str << std::endl;
}

void plotLayeredEllipsoids() {

    Graph2D *plt;
    plt = Graph2D::Instance();
    Vector radii(3);
    Vector center(3);
    radii << 50, 50, 50;
    center << -353, 0, 0;
    Mat X = GetEllipsoidXFromRadii(radii, 500);
    RandomRotatePoints(X);
    TranslatePoints(X, center);
    bool errorFlag;
    Ellipsoid flipsoid = GetVerifiedEllipsoid(X, radii, errorFlag);

    Mat TestX = GetEllipsoidPointShell(flipsoid, 250, 0.1);
    Mat ScanTestX = MakeLayeredX(TestX, flipsoid, 14);

    plt->scatter3d(X, "green");
    plt->scatter3d(ScanTestX, "red");
    plt->show();
}

void testScanningBias() {
    // See if distribution of a/c, a/b, b/c stays the same.

    Graph2D *plt;
    plt = Graph2D::Instance();

    int M = 20000; // Number of ellipsoids
    int N = 50;  // Number of available points
    double dh = 0.2; // Initial shell accuracy (this should not matter due to correction)
    int sliceSize = 5; // number of sections in z-direction the ellipsoid is sliced.
    double errorVariance = 10;

    double minRadii = 25;
    double maxRadii = 50;
    double transMax = 105;

    std::vector<double> gtac;
    std::vector<double> gtab;
    std::vector<double> gtbc;
    std::vector<double> ac;
    std::vector<double> ab;
    std::vector<double> bc;
    std::vector<double> Da;
    std::vector<double> minAngleAvg;
    double avgAvol = 0;
    double avgBvol = 0;
    double avgAatoz = 0;
    double avgActoz = 0;
    double avgBatoz = 0;
    double avgBctoz = 0;

    bool errorFlag;
    double error;
    int Total = M;

    while (M > 0) {
        std::cout << Total - M + 1 << "/" << Total << std::endl;
        Vector radii(3);
        Vector center(3);
        // Get random parameters from test parameters
        for (int i = 0; i < 3; ++i) {
            radii(i) = Random::randU(Random::generator) * (maxRadii - minRadii) + minRadii;
            center(i) = Random::randU(Random::generator) * transMax;
        }

        // Make a verified ellipsoid
        Mat X = GetEllipsoidXFromRadii(radii);
        RandomRotatePoints(X);
        Ellipsoid gtFlipsoid = GetVerifiedEllipsoid(X, radii, errorFlag);

        if (errorFlag) {
            continue;
        }

        // Make test points randomly from shell, layer them, add error
        Mat TestX = GetEllipsoidPointShell(gtFlipsoid, N, dh);
        Mat ScanTestX = MakeLayeredX(TestX, gtFlipsoid, sliceSize);
        TranslatePoints(ScanTestX, center);
        AddNormalDistError(ScanTestX, gtFlipsoid, errorVariance);

        // Fit the erroneous points
        Ellipsoid flipsoid = EllipsoidLeastSquareFit(ScanTestX, error, errorFlag);

        // Sort radii
        if (errorFlag) {
            continue;
        }

        EllipsoidComparisonResult result = compareEllipsoids(gtFlipsoid, flipsoid, 0, 0);
        gtac.push_back(result.Aac);
        gtab.push_back(result.Aab);
        gtbc.push_back(result.Abc);
        ac.push_back(result.Bac);
        ab.push_back(result.Bab);
        bc.push_back(result.Bbc);
        Da.push_back(result.Da);
        minAngleAvg.push_back(result.minAngleAvg);
        avgAvol += result.Avol;
        avgBvol += result.Bvol;
        avgAatoz += result.Aatoz;
        avgActoz += result.Actoz;
        avgBatoz += result.Batoz;
        avgBctoz += result.Bctoz;

        --M;

    }

    avgAvol /= Total;
    avgBvol /= Total;
    avgAatoz /= Total;
    avgActoz /= Total;
    avgBatoz /= Total;
    avgBctoz /= Total;

    plt->subplot(1, 2, 1);
    plt->title("Ground truth ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio a/c");
    plt->ylabel("frequency");
    plt->hist(gtac, "[i*0.125 + 1 for i in range(16)]", true);
    plt->subplot(1, 2, 2);
    plt->title("Fitted ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio a/c");
    plt->ylabel("frequency");
    plt->hist(ac, "[i*0.125 + 1 for i in range(16)]", true);
    plt->show();

    plt->subplot(1, 2, 1);
    plt->title("Ground truth ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio a/b");
    plt->ylabel("frequency");
    plt->hist(gtab, "[i*0.125 + 1 for i in range(16)]", true);
    plt->subplot(1, 2, 2);
    plt->title("Fitted ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio a/b");
    plt->ylabel("frequency");
    plt->hist(ab, "[i*0.125 + 1 for i in range(16)]", true);
    plt->show();

    plt->subplot(1, 2, 1);
    plt->title("Ground truth ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio b/c");
    plt->ylabel("frequency");
    plt->hist(gtbc, "[i*0.125 + 1 for i in range(16)]", true);
    plt->subplot(1, 2, 2);
    plt->title("Fitted ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio b/c");
    plt->ylabel("frequency");
    plt->hist(bc, "[i*0.125 + 1 for i in range(16)]", true);
    plt->show();

    plt->subplot(1, 2, 1);
    plt->title("Change in largest radii");
    //plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("change in a");
    plt->ylabel("frequency");
    plt->hist(Da, 20, true);//, "[i*0.125 + 1 for i in range(16)]", true);
    plt->subplot(1, 2, 2);
    plt->title("Average minimum radii angle to original ellipsoid");
    //plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("Avg min angle (deg)");
    plt->ylabel("frequency");
    plt->hist(minAngleAvg, 25, true);//, "[i*0.125 + 1 for i in range(16)]", true);
    plt->show();

    std::cout << "Average change in volume: " << avgBvol - avgAvol << " (" << (avgBvol - avgAvol)/avgAvol*100 << "%%)" << std::endl;
    std::cout << "Average min angle to z direction (deg) ref - a: " << avgAatoz << " c: " << avgActoz << std::endl;
    std::cout << "Average min angle to z direction (deg) fit - a: " << avgBatoz << " c: " << avgBctoz << std::endl;
}

void modelConsistency() {
    // See if distribution of a/c, a/b, b/c stays the same.

    Graph2D *plt;
    plt = Graph2D::Instance();

    int M = 20000; // Number of ellipsoids
    int N = 50;  // Number of available points
    double dh = 0.2; // Initial shell accuracy (this should not matter due to correction)
    double errorVariance = 10;

    double minRadii = 25;
    double maxRadii = 50;
    double transMax = 105;

    std::vector<double> gtac;
    std::vector<double> gtab;
    std::vector<double> gtbc;
    std::vector<double> ac;
    std::vector<double> ab;
    std::vector<double> bc;
    std::vector<double> Da;
    std::vector<double> minAngleAvg;
    double avgAvol = 0;
    double avgBvol = 0;

    bool errorFlag;
    double error;
    int Total = M;

    while (M > 0) {
        std::cout << Total - M + 1 << "/" << Total << std::endl;
        Vector radii(3);
        Vector center(3);
        // Get random parameters from test parameters
        for (int i = 0; i < 3; ++i) {
            radii(i) = Random::randU(Random::generator) * (maxRadii - minRadii) + minRadii;
            center(i) = Random::randU(Random::generator) * transMax;
        }

        // Make a verified ellipsoid
        Mat X = GetEllipsoidXFromRadii(radii);
        RandomRotatePoints(X);
        Ellipsoid gtFlipsoid = GetVerifiedEllipsoid(X, radii, errorFlag);

        if (errorFlag) {
            continue;
        }

        // Make test points randomly from shell, layer them, add error
        Mat TestX = GetEllipsoidPointShell(gtFlipsoid, N, dh);
        //Mat ScanTestX = MakeLayeredX(TestX, gtFlipsoid, 10, true);
        TranslatePoints(TestX, center);
        AddNormalDistError(TestX, gtFlipsoid, errorVariance, false);

        // Fit the erroneous points
        Ellipsoid flipsoid = EllipsoidLeastSquareFit(TestX, error, errorFlag);

        // Sort radii
        if (errorFlag) {
            continue;
        }

        EllipsoidComparisonResult result = compareEllipsoids(gtFlipsoid, flipsoid, 0, 0);
        gtac.push_back(result.Aac);
        gtab.push_back(result.Aab);
        gtbc.push_back(result.Abc);
        ac.push_back(result.Bac);
        ab.push_back(result.Bab);
        bc.push_back(result.Bbc);
        Da.push_back(result.Da);
        minAngleAvg.push_back(result.minAngleAvg);
        avgAvol += result.Avol;
        avgBvol += result.Bvol;

        --M;

    }

    avgAvol /= Total;
    avgBvol /= Total;

    plt->subplot(1, 2, 1);
    plt->title("Ground truth ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio a/c");
    plt->ylabel("frequency");
    plt->hist(gtac, "[i*0.125 + 1 for i in range(16)]", true);
    plt->subplot(1, 2, 2);
    plt->title("Fitted ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio a/c");
    plt->ylabel("frequency");
    plt->hist(ac, "[i*0.125 + 1 for i in range(16)]", true);
    plt->show();

    plt->subplot(1, 2, 1);
    plt->title("Ground truth ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio a/b");
    plt->ylabel("frequency");
    plt->hist(gtab, "[i*0.125 + 1 for i in range(16)]", true);
    plt->subplot(1, 2, 2);
    plt->title("Fitted ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio a/b");
    plt->ylabel("frequency");
    plt->hist(ab, "[i*0.125 + 1 for i in range(16)]", true);
    plt->show();

    plt->subplot(1, 2, 1);
    plt->title("Ground truth ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio b/c");
    plt->ylabel("frequency");
    plt->hist(gtbc, "[i*0.125 + 1 for i in range(16)]", true);
    plt->subplot(1, 2, 2);
    plt->title("Fitted ratio distribution");
    plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("ratio b/c");
    plt->ylabel("frequency");
    plt->hist(bc, "[i*0.125 + 1 for i in range(16)]", true);
    plt->show();

    plt->subplot(1, 2, 1);
    plt->title("Change in largest radii");
    //plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("change in a");
    plt->ylabel("frequency");
    plt->hist(Da, 20, true);//, "[i*0.125 + 1 for i in range(16)]", true);
    plt->subplot(1, 2, 2);
    plt->title("Average minimum radii angle to original ellipsoid");
    //plt->x_lim(1,3.25);
    plt->y_lim(0,0.5);
    plt->xlabel("Avg min angle (deg)");
    plt->ylabel("frequency");
    plt->hist(minAngleAvg, 25, true);//, "[i*0.125 + 1 for i in range(16)]", true);
    plt->show();

    std::cout << "Average change in volume: " << avgBvol - avgAvol << " (" << (avgBvol - avgAvol)/avgAvol*100 << "%%)" << std::endl;

}

void testAngleCorrelation() {

    Graph2D *plt;
    plt = Graph2D::Instance();

    int M = 1000000;

    std::vector<double> xes;
    std::vector<double> yes;
    std::vector<double> avs;

    xes = plt->linspace(10, 100, 9);
    avs = plt->linspace(0, PI/2/9*10, 10);
    for (int i = 0; i < 9; ++i) {
        //yes.push_back(std::cos(avs[i]) - std::cos(avs[i+1]));
        yes.push_back(std::sin(avs[i+1]) - std::sin(avs[i]));
    }

    Vector v(3);
    v << 1,0,0;
    Vector zUnit(3);
    zUnit << 0,0,1;

    std::vector<double> angles;

    for (int i = 0; i < M; ++i) {
        Mat R = RandomRotationMatrix();
        Vector x = R*v;
        double value = std::min(std::acos(x.dot(zUnit))/(2*PI)*360,std::acos(x.dot(-zUnit))/(2*PI)*360);
        angles.push_back(value);
    }
    plt->title("Minimum angle to z distribution");
    //plt->x_lim(0,90);
    plt->y_lim(0,0.20);
    plt->xlabel("Angle (deg)");
    plt->ylabel("frequency");
    plt->hist(angles, 9, true);
    plt->plot(xes, yes, "'red'", 2.5);
    plt->show();
}



void minimizerTest() {
    Minimizer minimizer;

    std::function<double(Vector const &)> F =
        [&](Vector const &v)
    {
        double a = 1;
        double b = 100;
        double x = v(0);
        double y = v(1);
        //return std::pow(x + 2*y - 7, 2) + std::pow(2*x + y - 5, 2);   //
        return std::pow(a - x,2) + b * std::pow(y - std::pow(x,2),2); // Rosenbrock x_min = (1,1)
        //return std::pow(x,2) + std::pow(y,2);                         // x^2 + y^2, x_min = (0,0)
    };

    Vector xInit(2);
    xInit << -1, 1;

    bool fail_flag = false;

    Vector xMin = minimizer.Minimize(F, xInit, fail_flag);
    std::cout << "xMin = " << xMin << std::endl;
}



void ellipsoidMinimizerTest() {

    EllipsoidMinimizer ellipMini;

    Vector trueCenter(3);
    Mat trueRotate(3,3);
    Vector trueRadii(3);

    trueRadii << 50, 35, 25;
    trueCenter << 3, 2, 1;

    Mat X = GetEllipsoidXFromRadii(trueRadii);
    trueRotate = RandomRotatePoints(X);

    TranslatePoints(X, trueCenter);

    Vector testCenter(3);
    Mat testRotate(3,3);
    Vector testRadii(3);

    testCenter << 1, 2, 3;
    testRotate = Mat::Identity(3,3);
    testRadii << 50, 50, 50;

    bool fail_flag = false;

    ellipMini.GetMinimum(X, testRadii, testCenter, testRotate, fail_flag);
    std::cout << "true radii = " << trueRadii.transpose() << std::endl;
    std::cout << "fitt radii = " << testRadii.transpose() << std::endl;
    std::cout << "true center = " << trueCenter.transpose() << std::endl;
    std::cout << "fitt center = " << testCenter.transpose() << std::endl;
    std::cout << "true rotation matrix:" << std::endl << trueRotate << std::endl;
    std::cout << "fitt rotation matrix:" << std::endl << testRotate << std::endl;

    // Test how to get major radii direction

}

void plotEnergyFootprint() {

    Graph2D *plt;
    plt = Graph2D::Instance();

    Vector trueCenter(3);
    Vector trueRadii(3);

    trueRadii << 50, 40, 30;
    trueCenter << 10, 20, -15;

    Mat X = GetEllipsoidXFromRadii(trueRadii);
    Mat R = RandomRotatePoints(X);

    for (int i = 0; i < 4; ++i) {
        plt->subplot(2, 2, i+1);
        plt->ylabel("Energy");
        Vector testCenter(3);
        Vector testRadii(3);
        Vector testQuat(4);
        testRadii << 25, 35, 50;
        testCenter << 10, 10, -15;
        testQuat = MatrixToQuat(R);
        std::string s = "Quat_2 (" + std::to_string(testQuat(1)) + ")";
        if (i == 0) {
            std::vector<double> xes = plt->linspace(-10,20,150);
            std::vector<double> yes;
            for (int i = 0; i < xes.size(); ++i) {
                testCenter(0) = xes[i];
                yes.push_back(EllipsoidEnergy(X, testRadii, testCenter, testQuat));
            }
            plt->xlabel("Center x (10)");
            plt->plot(xes, yes);
        } else if (i == 1) {
            std::vector<double> xes = plt->linspace(-20,65,150);
            std::vector<double> yes;
            for (int i = 0; i < xes.size(); ++i) {
                testRadii(0) = xes[i];
                yes.push_back(EllipsoidEnergy(X, testRadii, testCenter, testQuat));
            }
            plt->xlabel("Radii a (50)");
            plt->plot(xes, yes);
        }  else if (i == 2) {
            //std::vector<double> xes = plt->linspace(-10,10,1501);
            std::vector<double> xes = plt->linspace(0,1,1501);
            std::vector<double> yes;
            for (int i = 0; i < xes.size(); ++i) {
                testQuat(1) = xes[i];
                yes.push_back(EllipsoidEnergy(X, testRadii, testCenter, testQuat));
            }
            plt->xlabel(s);
            plt->plot(xes, yes);
        } else {
            std::vector<double> xes = plt->linspace(-50,50,201);
            std::vector<double> yes;
            for (int i = 0; i < xes.size(); ++i) {
                testCenter(0) = xes[i];
                testRadii(0) = xes[i];
                testQuat(1) = xes[i];
                yes.push_back(EllipsoidEnergy(X, testRadii, testCenter, testQuat));
            }
            plt->xlabel("Center x, Radii a, Quat_2");
            plt->plot(xes, yes);
        }
    }
    plt->show();
    //std::vector<double> xes;

    Vector testQuat = MatrixToQuat(R);
    std::cout << "Insanity check" << std::endl;
    std::cout << R*R.inverse() << std::endl;
    std::cout << "Real check" << std::endl;
    std::cout << QuatToMatrix(testQuat)*R.inverse() << std::endl;
}

void fitMahdiehEllipsoids() {

    Graph2D *plt;
    plt = Graph2D::Instance();

    std::vector<double> ED_averageResiduleDivMajorRadii;
    std::vector<double> ED_ac;
    std::vector<double> ED_area;
    std::vector<double> AD_averageResiduleDivMajorRadii;
    std::vector<double> AD_ac;
    std::vector<double> AD_area;

    TaskMaster<EllipsoidComparisonResult, std::string> taskMaster;

    std::function<EllipsoidComparisonResult(std::string const &)> FA =
        [&](std::string const &s)
    {
        Mat X = loadVesicleFromFolder(s);
        EllipsoidMinimizer ellipMini;
        bool errorFlag;
        double hyperError;
        double e1, e2;
        //Ellipsoid flipsoidA = ellipMini.FitEllipsoid(X, hyperError, errorFlag); // euclidean distance ellipsoid fit
        Ellipsoid flipsoidA = ellipMini.FitSphere(X, hyperError, errorFlag); // euclidean distance ellipsoid fit
        e1 = ellipMini.getAverageDist();
        Ellipsoid flipsoidB = ellipMini.AlgebraicDistanceEllipsoidFit(X, hyperError, errorFlag); // algebraic distance ellipsoid fit
        e2 = ellipMini.getAverageDist();
        EllipsoidComparisonResult cResult = compareEllipsoids(flipsoidA, flipsoidB, e1, e2);
        return cResult;
    };

    std::vector<std::string> params;

    for (int i = 1; i <= 86; ++i) {
        std::string folderNumber = std::to_string(i);
        std::string folderRootPath = "/home/dith/Dropbox/Studie Personlig/cpp_numerics/src/projects/EllipsoidFitting/data/Vesicle" + folderNumber;
        params.push_back(folderRootPath);
    }

    taskMaster.runJobs(FA, params);
    std::vector<EllipsoidComparisonResult> results = taskMaster.getResults();

    for (int i = 0; i < 86; ++i) {
        EllipsoidComparisonResult cResult = results[i];

        Vector radiiA = cResult.flipsoidA.getRadii();
        ED_averageResiduleDivMajorRadii.push_back(cResult.AAverageError/radiiA(0));
        ED_ac.push_back(cResult.Aac);
        ED_area.push_back(cResult.ASurfaceArea);

        Vector radiiB = cResult.flipsoidB.getRadii();
        AD_averageResiduleDivMajorRadii.push_back(cResult.BAverageError/radiiB(0));
        AD_ac.push_back(cResult.Bac);
        AD_area.push_back(cResult.BSurfaceArea);
    }

    // PLOTTING
    plt->subplot(1, 2, 1);
    plt->title("Average residule (major radii normalized)");
    //plt->x_lim(1,3.25);
    plt->y_lim(0,0.6);
    plt->xlabel("Normalized average residule");
    plt->ylabel("frequency");
    plt->hist(AD_averageResiduleDivMajorRadii, "[(i+1)*0.025 for i in range(20)]", true);
    plt->hist(ED_averageResiduleDivMajorRadii, "[(i+1)*0.025 for i in range(20)]", true, 0.7);
    plt->legend(std::vector<std::string>({"Algebraic","Euclidean"}), 1);
    plt->subplot(1, 2, 2);
    plt->title("Fitted ratio distribution");
    plt->x_lim(1,2.6);
    plt->y_lim(0,0.6);
    plt->xlabel("ratio a/c");
    plt->ylabel("frequency");
    plt->hist(AD_ac, "[i*0.125 + 1 - 0.125/2 for i in range(16)]", true);
    plt->hist(ED_ac, "[i*0.125 + 1 - 0.125/2 for i in range(16)]", true, 0.7);
    plt->legend(std::vector<std::string>({"Algebraic","Euclidean"}), 1);
    plt->show();


    plt->subplot(1, 2, 1);
    plt->title("Surface area distribution");
    //plt->x_lim(1,3.25);
    plt->y_lim(0,0.6);
    plt->xlabel("Surface area");
    plt->ylabel("frequency");
    plt->hist(AD_area, 12, true);
    plt->hist(ED_area, 12, true, 0.7);
    plt->subplot(1, 2, 2);
    plt->title("Surface area distribution");
    //plt->x_lim(1,3.25);
    plt->y_lim(0,0.6);
    plt->xlabel("Surface area");
    plt->ylabel("frequency");
    plt->hist(AD_area, 12, true);
    plt->hist(ED_area, 12, true, 0.7);
    plt->show();

}

void sphereNormalizationPlot() {
    Graph2D *plt;
    plt = Graph2D::Instance();

    plt->plotSphere3d();

    for (int i = 1; i <= 86; ++i) {

        //if (i != 5) {
        //    continue;
        //}
        // Init
        EllipsoidMinimizer ellipMini;
        bool errorFlag;
        double hyperError;

        // Load vesicle points
        std::string folderNumber = std::to_string(i);
        std::string folderRootPath = "/home/dith/Dropbox/Studie Personlig/cpp_numerics/src/projects/EllipsoidFitting/data/Vesicle" + folderNumber;

        Mat X = loadVesicleFromFolder(folderRootPath);

        Ellipsoid flipsoid = ellipMini.FitSphere(X, hyperError, errorFlag); // euclidean distance ellipsoid fit

        Vector radii = flipsoid.getRadii();
        TranslatePoints(X, -flipsoid.getCenter());
        X = X * 1.0/radii(0);

        plt->scatter3d(X, "green");
        //plt->plotTriSurfVesicle(X);


    }

    plt->show();
}

void rotationMatrixTest() {
    Vector mu(3);
    Vector radii(3);
    //Vector eulerAngles(3);
    //eulerAngles << PI/2.0, 0.2, 0.0;
    radii << 1.0/pow(25,2), 1.0/pow(35,2), 1.0/pow(50,2);
    mu << 0, 0, 0;
    Mat D = radii.asDiagonal();
    //Mat R = GetEulerRotationMatrix(eulerAngles);
    Mat R = RandomRotationMatrix();
    //R.transposeInPlace();
    //D.col(0).swap(D.col(1));
    //D.row(0).swap(D.row(1));
    //R.col(0).swap(R.col(2));
    //R.transposeInPlace();
    Mat P(3,3);
    P << 0, 1, 0,
         1, 0, 0,
         0, 0, 1;
    std::cout << "Diagonal Matrix:" << std::endl <<  D << std::endl;
    std::cout << "Rotation Matrix:" << std::endl <<  R << std::endl;
    Vector x(3);
    x << 25, 0, 0;
    x = R.transpose()*x;
    std::cout << "Permuted Diagonal Matrix:" << std::endl <<  P*D*P << std::endl;
    std::cout << "Permuted Rotation Matrix:" << std::endl <<  P*R << std::endl;

    std::cout << (x-mu).transpose()*(P*R).transpose()*P*D*P*P*R*(x-mu) - 1 << std::endl;

    Vector y = R.row(0)*25;
    R.row(2) = -R.row(2);
    std::cout << R << std::endl;
    std::cout << (y-mu).transpose()*R.transpose()*D*R*(y-mu) - 1 << std::endl;
}

void radii2DErrorPlotBothMethods(bool translate, bool rotate)
{
    double MajorRadii = 50;
    int numSimulations = 10000;

    int errorSteps = 100;
    double minError = 0;
    double maxError = 10;
    double derror = (maxError - minError) / errorSteps;


    std::vector<double> residualErrors(errorSteps);
    std::vector<double> radiiAError(errorSteps);
    std::vector<double> radiiBError(errorSteps);
    std::vector<double> radiiCError(errorSteps);
    std::vector<double> radiiAErrorVar(errorSteps);
    std::vector<double> radiiBErrorVar(errorSteps);
    std::vector<double> radiiCErrorVar(errorSteps);
    std::vector<double> acErrors(errorSteps);

    #pragma omp parallel for
    for (int m = 0; m < errorSteps; ++m)
    {
        std::vector<double> dataPrintA;
        std::vector<double> dataPrintB;
        std::vector<double> dataPrintC;

        double residualError = m * derror;
        double acError = 0;

        int n = 0;

        while (n < numSimulations)
        {
            EllipsoidMinimizer ellipMini;

            Vector trueCenter(3);
            Mat trueRotate(3,3);
            Vector trueRadii(3);

            trueRadii << MajorRadii
                       , Random::randU(Random::generator) * 20 + 30
                       , Random::randU(Random::generator) * 20 + 30;

            if (trueRadii(1) < trueRadii(2)) std::swap(trueRadii(1), trueRadii(2));

            if (rotate)
            {
                trueCenter << Random::randU(Random::generator)*100
                            , Random::randU(Random::generator)*100
                            , Random::randU(Random::generator)*100;
            }
            else
            {
                trueCenter << 0, 0, 0;
            }

            Mat X = GetEllipsoidXFromRadii(trueRadii);

            trueRotate = Mat::Identity(3,3);

            if (rotate) trueRotate = RandomRotatePoints(X);
            if (translate) TranslatePoints(X, trueCenter);

            Vector testCenter(3);
            Mat testRotate(3,3);
            Vector testRadii(3);

            testCenter << 0, 0, 0;
            testRotate = Mat::Identity(3,3);
            testRadii << 50, 50, 50;


            //std::cout << "true radii = " << trueRadii.transpose() << std::endl;
            //std::cout << "true center = " << trueCenter.transpose() << std::endl;
            //std::cout << "true rotation matrix:" << std::endl << trueRotate << std::endl;

            /*ellipMini.GetMinimum(X, testRadii, testCenter, testRotate);
            std::cout << "fitt radii = " << testRadii.transpose() << std::endl;
            std::cout << "fitt center = " << testCenter.transpose() << std::endl;
            std::cout << "fitt rotation matrix:" << std::endl << testRotate << std::endl;*/

            double hyperError;
            bool errorFlag;

            Ellipsoid gtFlipsoid = GetVerifiedEllipsoid(X, trueRadii, errorFlag);

            if (errorFlag) continue;

            AddNormalDistError(X, gtFlipsoid, residualError, false);

            Ellipsoid flipsoid = ellipMini.AlgebraicDistanceEllipsoidFit(X, hyperError, errorFlag); // algebraic distance ellipsoid fit
            if (errorFlag) continue;
            //std::cout << "fitt radii = " << flipsoid.getRadii().transpose() << std::endl;
            //std::cout << "fitt center = " << flipsoid.getCenter().transpose() << std::endl;
            //std::cout << "fitt rotation matrix:" << std::endl << flipsoid.getRotationMatrix() << std::endl;

            Vector radiiError = trueRadii - flipsoid.getRadii();
            for (int i = 0; i < 3; ++i)
            {
                //radiiError(i) = std::abs(radiiError(i));
            }

            dataPrintA.push_back(radiiError(0));
            dataPrintB.push_back(radiiError(1));
            dataPrintC.push_back(radiiError(2));
            acError += flipsoid.getRadii()(0)/flipsoid.getRadii()(2) - trueRadii(0)/trueRadii(2);

            ++n;
        }

        double EA = 0, EB = 0, EC = 0;
        double VA = 0, VB = 0, VC = 0;

        double N = dataPrintA.size();

        for (int i = 0; i < dataPrintA.size(); ++i)
        {
            EA += dataPrintA[i];
            EB += dataPrintB[i];
            EC += dataPrintC[i];
        }

        EA /= float(N);
        EB /= float(N);
        EC /= float(N);

        for (int i = 0; i < dataPrintA.size(); ++i)
        {
            VA += std::pow(dataPrintA[i] - EA, 2);
            VB += std::pow(dataPrintB[i] - EB, 2);
            VC += std::pow(dataPrintC[i] - EC, 2);
        }

        VA /= float(N-1);
        VB /= float(N-1);
        VC /= float(N-1);

        residualErrors[m] = residualError;

        radiiAError[m] = EA;
        radiiBError[m] = EB;
        radiiCError[m] = EC;

        radiiAErrorVar[m] = VA;
        radiiBErrorVar[m] = VB;
        radiiCErrorVar[m] = VC;

        acErrors[m] = acError / float(N);
    }

    std::cout << "res = [";
    for (int i = 0; i < errorSteps; ++i)
    {
        std::cout << residualErrors[i];
        if (i < errorSteps - 1) std::cout << ", ";
        else std::cout << "]" << std::endl;
    }

    std::cout << "ea = [";
    for (int i = 0; i < errorSteps; ++i)
    {
        std::cout << radiiAError[i];
        if (i < errorSteps - 1) std::cout << ", ";
        else std::cout << "]" << std::endl;
    }

    std::cout << "eb = [";
    for (int i = 0; i < errorSteps; ++i)
    {
        std::cout << radiiBError[i];
        if (i < errorSteps - 1) std::cout << ", ";
        else std::cout << "]" << std::endl;
    }

    std::cout << "ec = [";
    for (int i = 0; i < errorSteps; ++i)
    {
        std::cout << radiiCError[i];
        if (i < errorSteps - 1) std::cout << ", ";
        else std::cout << "]" << std::endl;
    }

    std::cout << "va = [";
    for (int i = 0; i < errorSteps; ++i)
    {
        std::cout << radiiAErrorVar[i];
        if (i < errorSteps - 1) std::cout << ", ";
        else std::cout << "]" << std::endl;
    }

    std::cout << "vb = [";
    for (int i = 0; i < errorSteps; ++i)
    {
        std::cout << radiiBErrorVar[i];
        if (i < errorSteps - 1) std::cout << ", ";
        else std::cout << "]" << std::endl;
    }

    std::cout << "vc = [";
    for (int i = 0; i < errorSteps; ++i)
    {
        std::cout << radiiCErrorVar[i];
        if (i < errorSteps - 1) std::cout << ", ";
        else std::cout << "]" << std::endl;
    }

    std::cout << "eac = [";
    for (int i = 0; i < errorSteps; ++i)
    {
        std::cout << acErrors[i];
        if (i < errorSteps - 1) std::cout << ", ";
        else std::cout << "]" << std::endl;
    }
}

void drift_effect_experiments_simulate_theoretic()
{
    double MajorRadii = 50;
    int numSimulations = 2500;//2500;

   /* std::vector<double> radii_kx1(numSimulations);
    std::vector<double> radii_ky1(numSimulations);
    std::vector<double> radii_kx2(numSimulations);
    std::vector<double> radii_ky2(numSimulations);*/

    std::vector<double> A(numSimulations);
    std::vector<double> B(numSimulations);
    std::vector<double> C(numSimulations);
    std::vector<double> D(numSimulations);
    std::vector<double> E(numSimulations);
    std::vector<double> F(numSimulations);

    //std::vector<double> Rj1(numSimulations);
    //std::vector<double> Rj2(numSimulations);
    //std::vector<double> Rj3(numSimulations);


    #pragma omp parallel for
    for (int n = 0; n < numSimulations; ++n)
    {
        Vector radii(3);
        Mat R(3,3);
        Mat H(3,3);
        Mat Di = Mat::Zero(3,3);

        radii << 2//5//MajorRadii
               , 1//3//Random::randU(Random::generator) * 20 + 30
               , 1;//2;//Random::randU(Random::generator) * 20 + 30;

        if (radii(1) < radii(2)) std::swap(radii(1), radii(2));

        for (int i = 0; i < 3; ++i)
        {
            Di(i,i) = 1.0/std::pow(radii(i), 2);
        }

        R = RandomRotationMatrix();

        H = R.transpose()*Di*R;

        //radii_kx1[n] = kx1;
        //radii_ky1[n] = ky1;
        //radii_kx2[n] = kx2;
        //radii_ky2[n] = ky2;

        A[n] = H(0,0);
        B[n] = H(1,1);
        C[n] = H(2,2);
        D[n] = H(0,1);
        E[n] = H(0,2);
        F[n] = H(1,2);

        //Rj1[n] = R(2,0);
        //Rj2[n] = R(2,1);
        //Rj3[n] = R(2,2);
    }

    /*std::cout << "kxs1 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << radii_kx1[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "kys1 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << radii_ky1[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "kxs2 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << radii_kx2[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "kys2 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << radii_ky2[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }*/

    std::cout << "A = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << A[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "B = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << B[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "C = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << C[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "D = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << D[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "E = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << E[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "F = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << F[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    /*std::cout << "Rj1 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << Rj1[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "Rj2 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << Rj2[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "Rj3 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << Rj3[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }*/
}

void drift_effect_experiments_simulate_with_fitting()
{
    double MajorRadii = 25;
    int numSimulations = 100;

    EllipsoidMinimizer ellipMini;

    std::vector<double> radii_kx1(numSimulations);
    std::vector<double> radii_ky1(numSimulations);
    std::vector<double> radii_kx2(numSimulations);
    std::vector<double> radii_ky2(numSimulations);

    std::mutex countMutex;

    int counter = 0;

    #pragma omp parallel for
    for (int n = 0; n < numSimulations; ++n)
    {
        Vector trueCenter(3);

        trueCenter << 0, 0, 0;

        Vector trueRadii(3);

        trueRadii << Random::randU(Random::generator) * 10 + 20
                   , Random::randU(Random::generator) * 10 + 20
                   , Random::randU(Random::generator) * 10 + 20;

        if (trueRadii(1) < trueRadii(2)) std::swap(trueRadii(1), trueRadii(2));
        if (trueRadii(2) < trueRadii(3)) std::swap(trueRadii(2), trueRadii(3));
        if (trueRadii(1) < trueRadii(2)) std::swap(trueRadii(1), trueRadii(2));

        Mat X = GetEllipsoidXFromRadii(trueRadii);

        Mat R(3,3);

        R = RandomRotatePoints(X);

        //AddDrift(X, 0.1, 0.0);
        //AddDrift(X, 0.05, -0.2);

        double hyperError;

        Vector testCenter(3);
        Vector testRadii(3);
        Mat testRotate(3,3);

        testCenter << 0, 0, 0;
        testRotate = Mat::Identity(3,3);
        testRadii << 25, 25, 25;

        //bool fail_flag = false;

        //ellipMini.GetMinimum(X, testRadii, testCenter, testRotate, fail_flag);

        //Ellipsoid ellipsoid = Ellipsoid(testCenter, testRadii, testRotate);

        double error;
        bool errorFlag = false;

        Ellipsoid ellipsoid = ellipMini.FitEllipsoid(X, error, errorFlag);

        double kx1, ky1, kx2, ky2;

        ellipsoid.getRadii_k(kx1, ky1, kx2, ky2);

        countMutex.lock();
        ++counter;
        std::cout << "Completed: " << counter << std::endl;
        countMutex.unlock();

        radii_kx1[n] = kx1;
        radii_ky1[n] = ky1;
        radii_kx2[n] = kx2;
        radii_ky2[n] = ky2;
    }

    std::cout << "kxs1 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << radii_kx1[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "kys1 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << radii_ky1[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "kxs2 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << radii_kx2[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "kys2 = [";
    for (int i = 0; i < numSimulations; ++i)
    {
        std::cout << radii_ky2[i];
        if (i == numSimulations - 1) std::cout << "]\n";
        else std::cout << ", ";
    }
}

void driftEstimationMahdieh()
{
    EllipsoidMinimizer ellipMini;

    std::vector<double> kxs(86);
    std::vector<double> kys(86);
    std::vector<double> E_z_begin(86);
    std::vector<double> E_z_end(86);

    #pragma omp parallel for
    for (int i = 1; i <= 86; ++i)
    {
        std::string folderNumber = std::to_string(i);
        std::string folderRootPath = "/home/dith/Dropbox/Studie Personlig/cpp_numerics/src/projects/EllipsoidFitting/data/Vesicle" + folderNumber;

        Mat X = loadVesicleFromFolder(folderRootPath);

        double z_begin = X(0,2);
        double z_end = X(0,2);

        for (int i = 1; i < X.rows(); ++i)
        {
            z_begin = std::min(z_begin, X(i,2));
            z_end = std::max(z_end, X(i,2));
        }

        //CenterAtOrigin(X); / don't do this, already done by ellipsoid minimizer

        //AddDrift(X, 0.1, 0.0);
        AddDrift(X, 0.05, -0.2);

        Vector fittedCenter(3);
        Vector fittedRadii(3);
        Mat fittedRotate(3,3);


        fittedCenter << 0, 0, 0;
        fittedRotate = Mat::Identity(3,3);
        fittedRadii << 25, 25, 25;

        //ellipMini.GetMinimum(X, fittedRadii, fittedCenter, fittedRotate);
        //Ellipsoid ellipsoid = Ellipsoid(fittedCenter, fittedRadii, fittedRotate);

        double error;
        bool errorFlag = false;

        Ellipsoid ellipsoid = ellipMini.FitEllipsoid(X, error, errorFlag);

        double kd,kx, ky;

        ellipsoid.getRadii_k(kd,kd,kx,ky);

        kxs[i-1] = kx;
        kys[i-1] = ky;

        E_z_begin[i-1] = z_begin;
        E_z_end[i-1] = z_end;
    }

    std::cout << "kxs = [";
    for (int i = 0; i < 86; ++i)
    {
        std::cout << kxs[i];
        if (i == 86 - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "kys = [";
    for (int i = 0; i < 86; ++i)
    {
        std::cout << kys[i];
        if (i == 86 - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "E_z_begin = [";
    for (int i = 0; i < 86; ++i)
    {
        std::cout << E_z_begin[i];
        if (i == 86 - 1) std::cout << "]\n";
        else std::cout << ", ";
    }

    std::cout << "E_z_end = [";
    for (int i = 0; i < 86; ++i)
    {
        std::cout << E_z_end[i];
        if (i == 86 - 1) std::cout << "]\n";
        else std::cout << ", ";
    }
}

int main() {

    //testAngleCorrelation();
    //testScanningBias();
    //modelConsistency();
    //plotEnergyFootprint();
    //ellipsoidMinimizerTest();
    //fitMahdiehEllipsoids();
    //sphereNormalizationPlot();
    //minimizerTest();
    //rotationMatrixTest();
    //radii2DErrorPlotBothMethods(true, true);
    //drift_effect_experiments_simulate_theoretic(); // begone pest

    //drift_effect_experiments_simulate_with_fitting();
    //driftEstimationMahdieh();

    return 0;
}