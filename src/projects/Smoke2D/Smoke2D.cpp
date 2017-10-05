// Solves the incompressible Navier Stokes equations in 2D
// Copyright 2012, Kenny Erleben, DIKU.
// Converted to C++ 2016, Hans J. T. Stephensen.


#include "Smoke2D.hpp"
#include <iostream>

Params::Params(double width, double height, int Iin, int Jin) : G(1, Iin + 2, 1, Jin + 2, Iin + 2, Jin + 2) {
    fps          = 60;                // Frames per second
    T            = 3;                 // Total simulated time
    dt           = 0.01;              // Time step size
    viscosity    = 0.0000001;         // Viscosity coefficient
    rel_tol      = 0.01;              // Relative tolerance for linear system solver
    abs_tol      = 0.01;              // Absolute tolerance for linear system solver
    max_iter     = 250;               // Maximum number of iterations of linear system solver
    drop         = 0.000625;          // Buoyancy density drop coefficient
    lift         = 0.25;              // Buoyancy temperature lift coefficient
    dx           = width / (Iin-1);   // Grid spacing between two grid nodes.
    dy           = height / (Jin-1);  // Grid spacing between two grid nodes.
    I            = Iin + 2;           // Add ghost nodes to domain.
    J            = Jin + 2;           // Add ghost nodes to domain.
    Total        = I * J;             // Complete size
    smoke_rate   = 1000;              // The rate of smoke added to the simulation
    vorticity    = 75;                // The amount of small scale that should be added by the vorticity forces.
    ScreenWidth  = 540;
    ScreenHeight = 540;
    eps          = 10e-10;
    Jacobi       = false;

    smoke_source_i = std::make_pair(I/2 - I/10, I/2 + I/10);
    smoke_source_j = std::make_pair(J/2 - J/10, J/2 + J/10);
}

Smoke2D::Smoke2D(double width, double height, int I, int J) : params(width, height, I, J) {
    tk.addVar("Buoyancy ");
    tk.addVar("Vorticity");
    tk.addVar("Diffusion");
    tk.addVar("Pres Proj");
    tk.addVar("Adv u & v");
    tk.addVar("Add smoke");
    tk.addVar("Adv smoke");

    // Construct sparse matrices and Cholesky factorization
    if (!params.Jacobi) {
        double beta1 = 1.0 / (params.dx * params.dx);
        double gamma1 = 1.0 / (params.dy * params.dy);
        double alpha1 = -2.0 * (beta1 + gamma1);

        std::vector<double> values1 {alpha1, beta1, gamma1};
        std::vector<int> indices1 {0, 1, params.J};

        SpMat K1 = matHelp.BandMatrixSparse(params.Total, values1, indices1);

        chol1 = new Eigen::SimplicialCholesky<SpMat>(K1);

        double gamma2 = - (params.dt * params.viscosity) / (params.dy * params.dy);
        double beta2  = - (params.dt * params.viscosity) / (params.dx * params.dx);
        double alpha2 = 1 - 2 * (beta2 + gamma2);

        std::vector<double> values2 {alpha2, beta2, gamma2};
        std::vector<int> indices2 {0, 1, params.J};

        SpMat K2 = matHelp.BandMatrixSparse(params.Total, values2, indices2);

        chol2 = new Eigen::SimplicialCholesky<SpMat>(K2);
    }
}

Smoke2D::~Smoke2D() {
    if (!params.Jacobi) {
        delete chol1;
        delete chol2;
    }
}

void Smoke2D::run() {

    Mat u = Mat::Zero(params.I, params.J);
    Mat v = Mat::Zero(params.I, params.J);
    Mat smoke = Mat::Zero(params.I, params.J);

    double T_wanted = params.T;

    while (T_wanted > 0) {

        std::cout << "Time left: " << T_wanted << std::endl;

        double dt_wanted = 1.0f / params.fps;

        double min_dist = std::min(params.dx, params.dy);

        double grid_speed = min_dist / dt_wanted;

        while (dt_wanted > 0) {

            //--- Apply a crude CFL condition as a rule of thumb to avoid
            //--- stability issues. For details read:
            //
            // @book{osher2003level,
            //  title={Level set methods and dynamic implicit surfaces},
            //  author={Osher, S. and Fedkiw, R.P.},
            //  isbn={9780387954820},
            //  lccn={70022436},
            //  series={Applied mathematical sciences},
            //  url={http://books.google.dk/books?id=SQQI2vqWR7gC},
            //  year={2003},
            //  publisher={Springer}
            //  }
            //
            //  Around equation 3.10 and 3.11
            //

            double speed = get_max_speed(u, v);
            double max_speed = std::max(speed, grid_speed);
            double cfl_dt = std::max(std::min(dt_wanted, 0.9 * (min_dist / max_speed)),10e-4);

            params.dt = cfl_dt;

            compute_fractional_step(u, v, smoke, params);

            dt_wanted -= cfl_dt;
        }

#ifdef USE_SDL
        all_smoke.push_back(smoke);
#endif
        T_wanted -= 1.0f / params.fps;
    }
    std::cout << "Time left: " << T_wanted << std::endl;
    tk.printTimes();

#ifdef USE_SDL
    sdl_init(params); // If using SDL for animation
    sdl_animation(params);
    sdl_destroy();
#endif
}

void Smoke2D::compute_fractional_step(Mat &u, Mat &v, Mat &smoke, Params& params) {
    // COMPUTE_FRACTIONAL_STEP: Solves an Navier Stokes equation iteration
    // for velocity and smoke field.
    // INPUT:
    //    u  -- The current value of the x-component of the velocity field.
    //    v  -- The current value of the y-component of the velocity field.
    // smoke -- The current value of the smoke density field.
    // Output:
    //    u  -- The updated value of the x-component of the velocity field.
    //    v  -- The updated value of the y-component of the velocity field.
    // smoke -- The updated value of the smoke density field.
    //

    //--- First phase we time integrate the velocity field of the fluid ------

    bool timeit = true;



    Mat fv = Mat::Zero(params.I, params.J);
    Mat fu = Mat::Zero(params.I, params.J);

    //--- First step: integrate buoyancy force -------------------------------
    double timeStart;
    double timeEnd;
    if (timeit) timeStart = getCPUTime();
    compute_buoyancy_force(fv, smoke, params);
    if (timeit) timeEnd = getCPUTime();
    if (timeit) tk.addTime(0, timeEnd - timeStart);

    v += fv * params.dt;
    //--- Second step integrate vorticity force ------------------------------
    if (timeit) timeStart = getCPUTime();
    compute_vorticity_confinement(fu, fv, u, v, params);
    if (timeit) timeEnd = getCPUTime();
    if (timeit) tk.addTime(1, timeEnd - timeStart);

    v += fv * params.dt;
    u += fu * params.dt;

    //--- Third step: Integrate diffusion term -------------------------------

    Mat u0 = u;
    Mat v0 = v;
    if (timeit) timeStart = getCPUTime();
    compute_diffusion(u, u0, params);
    compute_diffusion(v, v0, params);
    if (timeit) timeEnd = getCPUTime();
    if (timeit) tk.addTime(2, timeEnd - timeStart);

    //--- Intermediate Step: extra projection for better accuracy ------------

    if (timeit) timeStart = getCPUTime();
    compute_pressure_projection(u, v, params);
    if (timeit) timeEnd = getCPUTime();
    if (timeit) tk.addTime(3, timeEnd - timeStart);

    //--- Fourth Step: Integrate advection term ------------------------------
    u0 = u;
    v0 = v;

    if (timeit) timeStart = getCPUTime();
    compute_advection(1, u, u0, v0, params);
    compute_advection(2, v, u0, v0, params);
    if (timeit) timeEnd = getCPUTime();
    if (timeit) tk.addTime(4, timeEnd - timeStart);


    //--- Fifth Step: Finaly projection step to make fluid divergence free ---
    if (timeit) timeStart = getCPUTime();
    compute_pressure_projection(u, v, params);
    if (timeit) timeEnd = getCPUTime();
    if (timeit) tk.addTime(3, timeEnd - timeStart);

    //--- Second phase, we time integrate our trace (smoke) density field ----

    //--- First step: We add some more new smoke -----------------------------

    if (timeit) timeStart = getCPUTime();
    add_smoke_at_source(smoke, params);
    if (timeit) timeEnd = getCPUTime();
    if (timeit) tk.addTime(5, timeEnd - timeStart);

    //--- Second step: We advect the smoke by the fluid velocity -------------
    if (timeit) timeStart = getCPUTime();
    compute_advection(0, smoke, u, v, params);
    if (timeit) timeEnd = getCPUTime();
    if (timeit) tk.addTime(6, timeEnd - timeStart);
}

void Smoke2D::compute_buoyancy_force(Mat &fv, Mat &smoke, Params &params) {
    double T_ambient = smoke.sum() / (params.Total);
    Mat T_ambient_mat = T_ambient * Mat::Ones(params.I, params.J);
    fv = -params.drop * smoke
       +  params.lift * (smoke - T_ambient_mat);
}

void Smoke2D::compute_vorticity_confinement(Mat &fu, Mat &fv, Mat &u, Mat &v, Params &params) {

    //--- w =  nabla times (u,v)^T -------------------------------------------
    Mat w = Mat::Zero(params.I, params.J);
    int inner_I = params.I - 2;
    int inner_J = params.J - 2;

    w.block(1, 1, inner_I, inner_J) = (v.block(2, 1, inner_I, inner_J) - v.block(0, 1, inner_I, inner_J)) * 1.0 / (2.0 * params.dx)
                          - (u.block(1, 2, inner_I, inner_J) - u.block(1, 0, inner_I, inner_J)) * 1.0 / (2.0 * params.dy);

    set_boundary_conditions(0, w, params);

    //--- n = nabla |w| ------------------------------------------------------

    Mat nu = Mat::Zero(params.I, params.J);
    Mat nv = Mat::Zero(params.I, params.J);

    Mat aw(params.I, params.J);
    for (int i = 0; i < params.I; ++i) {
        for (int j = 0; j < params.J; ++j) {
            aw(j, i) = std::abs(w(j, i));
        }
    }

    nu.block(1, 1, inner_I, inner_J) = (aw.block(2, 1, inner_I, inner_J) - aw.block(0, 1, inner_I, inner_J)) * 1.0 / (2.0 * params.dx);
    nv.block(1, 1, inner_I, inner_J) = (aw.block(1, 2, inner_I, inner_J) - aw.block(1, 0, inner_I, inner_J)) * 1.0 / (2.0 * params.dy);

    set_boundary_conditions(0, nu, params);
    set_boundary_conditions(0, nv, params);

    //--- N = n / | n | ------------------------------------------------------
    Mat norm_n = nu.array() * nu.array() + nv.array() * nv.array();

    Mat divn = Mat::Zero(params.I, params.J);

    for (int i = 0; i < params.I; ++i) {
        for (int j = 0; j < params.J; ++j) {
            divn(j,i) = sqrt(norm_n(j,i)) + params.eps;
        }
    }

    nu = nu.array() / divn.array();
    nv = nv.array() / divn.array();

    //--- f = N x w ----------------------------------------------------------
    double h = std::min(params.dx, params.dy);
    fu = (nv.array() *  w.array()) * (params.vorticity * h);
    fv = (nu.array() * -w.array()) * (params.vorticity * h);
}

void Smoke2D::compute_diffusion(Mat &x, Mat x0, Params &params) {


    //if (params.Jacobi) {
        double gamma = - (params.dt * params.viscosity) / (params.dy * params.dy);
        double beta  = - (params.dt * params.viscosity) / (params.dx * params.dx);
        double alpha = 1 - 2 * (beta + gamma);
        jacobi_iterative_solver(0, alpha, beta, gamma, x, x0, params);
    /*} else {
        x = sparse_solver(2, x0);
        x.resize(params.I, params.J);
    }*/
}

void Smoke2D::compute_advection(int type, Mat &phi, Mat &u, Mat &v, Params &params) {

    Mat y_cur = params.G.getAllX();
    Mat x_cur = params.G.getAllY();

    Mat x_old = x_cur - u * (params.dt / params.dx);
    Mat y_old = y_cur - v * (params.dt / params.dx);

    for (int j = 0; j < params.J; ++j) {
        for (int i = 0; i < params.I; ++i) {
            if (x_cur(i, j) < 1) x_cur(i, j) = 1;
            if (x_cur(i, j) > params.I) x_cur(i, j) = params.I;
            if (y_cur(i, j) < 1) y_cur(i, j) = 1;
            if (y_cur(i, j) > params.J) y_cur(i, j) = params.J;
        }
    }

    x_old.resize(params.Total, 1);
    y_old.resize(params.Total, 1);

    Mat new_phi = matHelp.interpolate(params.G, phi, y_old, x_old);

    new_phi.resize(params.I, params.J);

    phi = new_phi;

    set_boundary_conditions(type, phi, params);
}

void Smoke2D::compute_pressure_projection(Mat &u, Mat &v, Params &params) {
    Mat b = Mat::Zero(params.I, params.J);
    int inner_I = params.I - 2;
    int inner_J = params.J - 2;

    b.block(1, 1, inner_I, inner_J) = (u.block(2, 1, inner_I, inner_J) - u.block(0, 1, inner_I, inner_J)) * 1.0 / (2.0 * params.dx)
                                    + (v.block(1, 2, inner_I, inner_J) - v.block(1, 0, inner_I, inner_J)) * 1.0 / (2.0 * params.dy);

    set_boundary_conditions(0, b, params);

    Mat p;

    if (params.Jacobi) {
        double beta = 1.0 / (params.dx * params.dx);
        double gamma = 1.0 / (params.dy * params.dy);
        double alpha = -2.0 * (beta + gamma);
        p = Mat::Zero(params.I, params.J);
        jacobi_iterative_solver(0, alpha, beta, gamma, p, b, params);
    } else {
        p = sparse_solver(1, b); // matAlg.SparseSolve(K, b, params.eps, params.max_iter, 2);
        p.resize(params.I, params.J);
    }

    Mat nabla_px = Mat::Zero(inner_I, inner_J);
    Mat nabla_py = Mat::Zero(inner_I, inner_J);

    nabla_px = (p.block(2, 1, inner_I, inner_J) - p.block(0, 1, inner_I, inner_J)) * 1.0 / (2.0 * params.dx);
    nabla_py = (p.block(1, 2, inner_I, inner_J) - p.block(1, 0, inner_I, inner_J)) * 1.0 / (2.0 * params.dy);

    u.block(1, 1, inner_I, inner_J) -= nabla_px;
    v.block(1, 1, inner_I, inner_J) -= nabla_py;

    set_boundary_conditions(1, u, params);
    set_boundary_conditions(2, v, params);
}

void Smoke2D::set_boundary_conditions(int type, Mat &phi, Params &params) {

    // SET_BOUNDARY_CONDITIONS: Sets the value of the ghost nodes such that
    // boundary conditions are satisfied.
    //
    // We always have a boxed domain where a static solid wall is placed between
    // the outer one ring ghost nodes and inner domain nodes.
    //
    // INPUT
    //
    //    type  -  Field type indicator, can be 0 (pressure field),
    //             1 (u-component of velocity field), 2(v-component
    //             of velocity field)
    //
    //    phi   - The field where boundaryc ocnditions should be applied to.
    //   params - The parameter values.
    //
    // OUTPUT
    //
    //      phi - The updated field wih ghost values set appropriately.
    //
    // Copyright 2012, Kenny Erleben, DIKU

    int I = params.I;
    int J = params.J;

    int inner_I = params.I - 2;
    int inner_J = params.J - 2;

    //------------------------------------------------------------------------
    // With out loss of generality let us take a look at a left solid wall. The
    // ghost node to the left is given by index i-1 and the fluid cell just
    // right of the wall is given by index i.
    //
    // A pressure value at the wall would be estimated using a CD scheme
    //
    //  p_wall =  (p_{i-1} + p_i)/2
    //
    //  D(p_wall,\vec n) = \nabla p_wall^T \vec n \approx   (p_i - p_{i-1})/ dx
    //
    // Similar formulas hold for u_wall and v_wall.

    if (type > 0) {

        //--- No slip condition means that velocity is zero at walls, (v_wall,u_wall)^T = 0
        phi.block(1, 0, inner_I, 1)     = - phi.block(1, 1, inner_I, 1);
        phi.block(1, J - 1, inner_I, 1) = - phi.block(1, J - 2, inner_I, 1);

        phi.block(0, 1, 1, inner_J)     = - phi.block(1, 1, 1, inner_J);
        phi.block(I - 1, 1, 1, inner_J) = - phi.block(I - 2, 1, 1, inner_J);

    } else {

        //--- Pressure difference across a solid wall is zero, dp_wall/dn = 0
        phi.block(1, 0, inner_I, 1)     = phi.block(1, 1, inner_I, 1);
        phi.block(1, J - 1, inner_I, 1) = phi.block(1, J - 2, inner_I, 1);

        phi.block(0, 1, 1, inner_J)     = phi.block(1, 1, 1, inner_J);
        phi.block(I - 1, 1, 1, inner_J) = phi.block(I - 2, 1, 1, inner_J);

    }

    //--- Finally we fix ghost corners ---------------------------------------
    phi(0,0)     = 0.5 * (phi(1, 0)     + phi(0, 1));
    phi(0,J-1)   = 0.5 * (phi(1, J-1)   + phi(0, J-2));
    phi(I-1,0)   = 0.5 * (phi(I-1, 1)   + phi(I-2, 0));
    phi(I-1,J-1) = 0.5 * (phi(I-2, J-1) + phi(I-1, J-2));
}

double Smoke2D::get_max_speed(Mat u, Mat v) {
    double result = 0.0f;

    for (int i = 0; i < u.cols(); ++i) {
        for (int j = 0; j < u.rows(); ++j) {
            double new_speed = sqrt( pow(u(j,i), 2) + pow(v(j,i), 2) );
            result = std::max(new_speed, result);
        }
    }

    return result;
}

void Smoke2D::add_smoke_at_source(Mat &smoke, Params &params) {

    for (int i = params.smoke_source_i.first; i < params.smoke_source_i.second; ++i) {
        for (int j = params.smoke_source_j.first; j < params.smoke_source_j.second; ++j) {
            smoke(i,j) += params.smoke_rate * params.dt;
        }
    }
}

void Smoke2D::jacobi_iterative_solver(int type, double alpha, double beta, double gamma, Mat &x, Mat b, Params &params) {

    int inner_I = params.I - 2;
    int inner_J = params.J - 2;

    for (int i = 0; i < params.max_iter; ++i) {

        set_boundary_conditions(type, x, params);

        Mat x_old = x;
        x.block(1, 1, inner_I, inner_J) = ( b.block(1, 1, inner_I, inner_J)
            -  beta * (x.block(2, 1, inner_I, inner_J) + x.block(0, 1, inner_I, inner_J))
            - gamma * (x.block(1, 2, inner_I, inner_J) + x.block(1, 0, inner_I, inner_J))
            ) * 1.0 / alpha;
        if ((x_old - x).lpNorm<Eigen::Infinity>() < params.eps) {
            break;
        }
    }

    set_boundary_conditions(type, x, params);
}

Vector Smoke2D::sparse_solver(int type, Vector b) {

    Vector x;

    if (type == 1) {
        x = chol1->solve(b);
    } else if (type == 2) {
        x = chol2->solve(b);
    }
    return x;
}

void Smoke2D::sdl_init(Params &params) {
#ifdef USE_SDL

    if (SDL_Init(SDL_INIT_VIDEO) != 0){
        std::cout << "E: DEBUG: SDL_Init ERROR: " << SDL_GetError() << std::endl;
        exit(1);
    } else {
        std::cout << "I: DEBUG: SDL_Init SUCCESS:" << std::endl;
    }

    win = SDL_CreateWindow("Smoke 2D", 100, 100, params.ScreenWidth, params.ScreenHeight, SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    if (win == nullptr){
        std::cout << "E: DEBUG: SDL_CreateWindow ERROR:" << SDL_GetError() << std::endl;
        SDL_Quit();
        exit(1);
    } else {
        std::cout << "I: DEBUG: SDL_CreateWindow SUCCESS:" << std::endl;
    }

    ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_SOFTWARE); //, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC
    if (ren == nullptr){
        SDL_DestroyWindow(win);
        std::cout << "E: DEBUG: SDL_CreateRenderer ERROR: " << SDL_GetError() << std::endl;
        SDL_Quit();
        exit(1);
    } else {
        std::cout << "I: DEBUG: SDL_CreateRenderer SUCCESS: " << std::endl;
    }
#endif
}

void Smoke2D::sdl_destroy() {
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
}

void Smoke2D::sdl_animation(Params &params) {
#ifdef USE_SDL

    int frame_index = 0;
    int N = all_smoke.size();
    bool quit = false;
    while (!quit) {

        SDL_Event event;
        if (SDL_PollEvent( &event )
            && event.type == SDL_KEYDOWN
            && event.key.keysym.sym == SDLK_ESCAPE) {
            quit = true;
        }

        Mat smoke = all_smoke[frame_index % N];

        SDL_RenderClear(ren);

        for (int j = 0; j < params.ScreenWidth; ++j) {
            for (int i = 0; i < params.ScreenHeight; ++i) {
                int color_val = smoke(i * params.I / params.ScreenWidth, j * params.J / params.ScreenHeight);
                color_val = std::min(color_val, 255);
                SDL_SetRenderDrawColor(ren, color_val, color_val, color_val, 255);
                SDL_RenderDrawPoint(ren, i, params.ScreenWidth - 1 - j); //params.ScreenHeight - 1 -
            }
        }

        SDL_SetRenderDrawColor(ren, 255, 0, 0, 0);
        SDL_RenderDrawLine(ren, 0, 0, params.ScreenWidth-1, 0);
        SDL_RenderDrawLine(ren, 0, 0, 0, params.ScreenHeight-1);
        SDL_RenderDrawLine(ren, params.ScreenWidth-1, 0, params.ScreenWidth-1, params.ScreenHeight-1);
        SDL_RenderDrawLine(ren, params.ScreenWidth-1, params.ScreenHeight-1, 0, params.ScreenHeight-1);

        SDL_SetRenderDrawColor(ren, 0, 0, 0, 0);

        SDL_RenderPresent(ren);

        // This is not the wanted fps - correct me
        SDL_Delay(80);
        ++frame_index;
    }
#endif
}
