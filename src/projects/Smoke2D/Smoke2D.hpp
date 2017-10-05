#ifndef SMOKE_2D_HPP
#define SMOKE_2D_HPP

#include "../../functions/types.hpp"
#include "../../functions/MatrixHelpers.hpp"
#include "../../functions/MatrixAlgorithms.hpp"
#include "../../debug/benchmark.hpp"

#ifdef USE_SDL
    #include <SDL.h>
#endif

#include <vector>

class Params {
public:

    Params(double width, double height, int I, int J);

    int fps;
    double T;
    double dt;
    double viscosity;
    double rel_tol;
    double abs_tol;
    int max_iter;
    double drop;
    double lift;
    double dx;
    double dy;
    int I;
    int J;
    int Total;
    double smoke_rate;
    double vorticity;
    int ScreenWidth;
    int ScreenHeight;
    double eps;
    bool Jacobi;

    std::pair<int, int> smoke_source_i;
    std::pair<int, int> smoke_source_j;

    Meshgrid G;
};

class Smoke2D {
public:
    Smoke2D(double width, double height, int I, int J);
    ~Smoke2D();
    void run();
private:

    MatrixHelpers matHelp;
    MatrixAlgorithms matAlg;
    Eigen::SimplicialCholesky<SpMat> *chol1;
    Eigen::SimplicialCholesky<SpMat> *chol2;

    void compute_fractional_step(Mat &u, Mat &v, Mat &smoke, Params& params);

    void compute_buoyancy_force(Mat &fv, Mat &smoke, Params &params);

    void compute_vorticity_confinement(Mat &fu, Mat &fv, Mat &u, Mat &v, Params &params);

    void compute_diffusion(Mat &x, Mat x0, Params &params);

    void compute_advection(int type, Mat &u, Mat &u0, Mat &v0, Params &params);

    void compute_pressure_projection(Mat &u, Mat &v, Params &params);

    void set_boundary_conditions(int type, Mat &phi, Params &params);

    double get_max_speed(Mat u, Mat v);

    void add_smoke_at_source(Mat &smoke, Params &params);

    void jacobi_iterative_solver(int type, double alpha, double beta, double gamma, Mat &x, Mat b, Params &params);

    Vector sparse_solver(int type, Vector b);

    void sdl_init(Params &params);

    void sdl_destroy();

    void sdl_animation(Params &params);

#ifdef USE_SDL
    SDL_Window *win;
    SDL_Renderer *ren;
    std::vector<Mat> all_smoke;
    KeeperOfTimes tk;
    Params params;

#endif
};



#endif // SMOKE_2D_HPP