#ifndef CHANVESE3D_CU_H
#define CHANVESE3D_CU_H

#include "types.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

#define gIndex(z,y,x,ny,nx,s) (z+s)*(ny+2*s)*(nx+2*s) + (y+s)*(nx+2*s) + x+s

class ChanVese3D {
public:

    ChanVese3D() {
        dx = 1.0;
        numericEps = 0.000000000001;
        maxCurvature = 1.0/dx;
    };

    void segmentImageCPU(const double * const image, double * const Phi,
                      int nz, int ny, int nx) {

        int N = nz * ny * nx;

        double mu = 2500.0;//0.1*255*255;
        double nu = 1.0;
        double eps = 1.0;
        double lam1 = 1.00;
        double lam2 = 1.00;
        int maxIter = 25000;

        initPhiSphere(Phi, nz, ny, nx, nz/2, ny/2, nx/2, 5.0);
        //reinitializePhi(Phi, nz, ny, nx);
        double *C = (double*) malloc(N * sizeof(double));
        double *V = (double*) malloc(N * sizeof(double));

        for (int n = 0; n < maxIter; ++n) {

            double c1 = areaInside(image, Phi, nz, ny, nx, eps);
            double c2 = areaOutside(image, Phi, nz, ny, nx, eps);

            curvature(Phi, C, nz, ny, nx);

            #pragma omp parallel for collapse(3)
            for (int z = 0; z < nz; ++z) {
                for (int y = 0; y < ny; ++y) {
                    for (int x = 0; x < nx; ++x) {

                        V[gIndex(z,y,x,ny,nx,0)] = delta(Phi[gIndex(z,y,x,ny,nx,1)], eps)
                            * ( mu * C[gIndex(z,y,x,ny,nx,0)] - nu
                            - lam1 * pow(image[gIndex(z,y,x,ny,nx,0)] - c1, 2)
                            + lam2 * pow(image[gIndex(z,y,x,ny,nx,0)] - c2, 2)
                            );
                    }
                }
            }

            double Vmax = 0;
            for (int i = 0; i < nx*ny*nz; ++i) {
                Vmax = std::max(Vmax, fabs(V[i]));
            }

            double dt = 0.5*dx/Vmax;

            #pragma omp parallel for collapse(3)
            for (int z = 0; z < nz; ++z) {
                for (int y = 0; y < ny; ++y) {
                    for (int x = 0; x < nx; ++x) {
                        Phi[gIndex(z,y,x,ny,nx,1)] = Phi[gIndex(z,y,x,ny,nx,1)] + dt * V[gIndex(z,y,x,ny,nx,0)];
                    }
                }
            }

            applyBoundaryCondition(Phi, nz, ny, nx);

            //if (n > 0 and n % 100 == 0 and n < 0.9*maxIter) reinitializePhi(Phi, nz, ny, nx, 5);
            //reinitializePhi(Phi, nz, ny, nx);

            if (n % 250 == 0 or n < 50) std::cout << "n: " << n << " dt: " << dt << " vMax: " << Vmax << " averages: " << c1 << ", " << c2 << std::endl;
        }
        //reinitializePhi(Phi, nz, ny, nx);
        free(C);

    };

    void initPhiBubble(double * const out, const int nz, const int ny, const int nx, const double sz, const double sy, const double sx, const double bubblePeriod) {
        #pragma omp parallel for collapse(3)
        for (int z = -1; z < nz + 1; ++z) {
            for (int y = -1; y < ny + 1; ++y) {
                for (int x = -1; x < nx + 1; ++x) {
                    if (z == -1 or y == -1 or x == -1 or z == nz or y == ny or x == nx) {
                        out[gIndex(z,y,x,ny,nx,1)] = 0;
                    } else {
                        out[gIndex(z,y,x,ny,nx,1)] = 5
                            * cos(PI*double(x - sx)/bubblePeriod)
                            * cos(PI*double(y - sy)/bubblePeriod)
                            * cos(PI*double(z - sz)/bubblePeriod);
                    }
                }
            }
        }
    }

    void initPhiSphere(double * const out, const int nz, const int ny, const int nx, const double sz, const double sy, const double sx, const double radius) {
        #pragma omp parallel for collapse(3)
        for (int z = -1; z < nz + 1; ++z) {
            for (int y = -1; y < ny + 1; ++y) {
                for (int x = -1; x < nx + 1; ++x) {
                    //if (z == -1 or y == -1 or x == -1 or z == nz or y == ny or x == nx) {
                    //    out[gIndex(z,y,x,ny,nx,1)] = 0;
                    //} else {
                        out[gIndex(z,y,x,ny,nx,1)] = pow(double(x - sx), 2.0)
                                                   + pow(double(y - sy), 2.0)
                                                   + pow(double(z - sz), 2.0)
                                                   - radius * radius;
                    //}
                }
            }
        }
    }

    void curvature(const double * const Phi, double * const C, const int nz, const int ny, const int nx) {

        #pragma omp parallel for collapse(3)
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {

                    double Dx  = ( Phi[gIndex(z,  y,  x+1,ny,nx,1)]
                              -   Phi[gIndex(z,  y,  x-1,ny,nx,1)])/(2*dx);
                    double Dy  = ( Phi[gIndex(z,  y+1,x,  ny,nx,1)]
                              -   Phi[gIndex(z,  y-1,x,  ny,nx,1)])/(2*dx);
                    double Dz  = ( Phi[gIndex(z+1,y,  x,  ny,nx,1)]
                              -   Phi[gIndex(z-1,y,  x,  ny,nx,1)])/(2*dx);

                    double Dxx = ( Phi[gIndex(z,  y,  x+1,ny,nx,1)]
                              - 2*Phi[gIndex(z,  y,  x,  ny,nx,1)]
                              +   Phi[gIndex(z,  y,  x-1,ny,nx,1)])/(dx*dx);
                    double Dyy = ( Phi[gIndex(z,  y+1,x,  ny,nx,1)]
                              - 2*Phi[gIndex(z,  y,  x,  ny,nx,1)]
                              +   Phi[gIndex(z,  y-1,x,  ny,nx,1)])/(dx*dx);
                    double Dzz = ( Phi[gIndex(z+1,y,  x,  ny,nx,1)]
                              - 2*Phi[gIndex(z,  y,  x,  ny,nx,1)]
                              +   Phi[gIndex(z-1,y,  x,  ny,nx,1)])/(dx*dx);

                    double Dxy = (Phi[gIndex(z,  y+1,x+1,ny,nx,1)]
                              -  Phi[gIndex(z,  y+1,x-1,ny,nx,1)]
                              -  Phi[gIndex(z,  y-1,x+1,ny,nx,1)]
                              +  Phi[gIndex(z,  y-1,x-1,ny,nx,1)])/(4*dx*dx);
                    double Dxz = (Phi[gIndex(z+1,y,  x+1,ny,nx,1)]
                              -  Phi[gIndex(z+1,y,  x-1,ny,nx,1)]
                              -  Phi[gIndex(z-1,y,  x+1,ny,nx,1)]
                              +  Phi[gIndex(z-1,y,  x-1,ny,nx,1)])/(4*dx*dx);
                    double Dyz = (Phi[gIndex(z+1,y+1,x,  ny,nx,1)]
                              -  Phi[gIndex(z+1,y-1,x,  ny,nx,1)]
                              -  Phi[gIndex(z-1,y+1,x,  ny,nx,1)]
                              +  Phi[gIndex(z-1,y-1,x,  ny,nx,1)])/(4*dx*dx);

                    double G = pow(sqrt(Dx*Dx + Dy*Dy + Dz*Dz), 3);

                    C[gIndex(z,y,x,ny,nx,0)] = ( Dx*Dx*(Dyy + Dzz)
                                               + Dy*Dy*(Dxx + Dzz)
                                               + Dz*Dz*(Dxx + Dyy)
                                               - 2*Dx*Dy*Dxy
                                               - 2*Dx*Dz*Dxz
                                               - 2*Dy*Dz*Dyz
                                               ) / (G + numericEps);

                    // Clamp curvature to prevent excessive curvature of Phi
                    C[gIndex(z,y,x,ny,nx,0)] =
                        std::max(-maxCurvature,
                            std::min(C[gIndex(z,y,x,ny,nx,0)],
                                maxCurvature
                                )
                            );
                }
            }
        }
    }

    void reinitializePhi(double * const Phi,
                         const int nz,
                         const int ny,
                         const int nx,
                         const int numIter) {

        int N = nz*ny*nx;

        double *Phi0 = (double*) malloc(N * sizeof(double));

        #pragma omp parallel for collapse(3)
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    Phi0[gIndex(z,y,x,ny,nx,0)] = Phi[gIndex(z,y,x,ny,nx,1)];
                }
            }
        }


        for (int n = 0; n < numIter; ++n) {

            double *V = (double*) malloc(N * sizeof(double));

            #pragma omp parallel for collapse(3)
            for (int z = 0; z < nz; ++z) {
                for (int y = 0; y < ny; ++y) {
                    for (int x = 0; x < nx; ++x) {

                        double Phixyz = Phi[gIndex(z,y,x,ny,nx,1)];
                        double Phi0xyz = Phi0[gIndex(z,y,x,ny,nx,0)];

                        double a = ( Phixyz
                                -   Phi[gIndex(z,  y,  x-1,ny,nx,1)])
                                /   dx;
                        double b = ( Phi[gIndex(z,  y,  x+1,ny,nx,1)]
                                -   Phixyz)
                                /   dx;
                        double c = ( Phixyz
                                -   Phi[gIndex(z,  y-1,x,  ny,nx,1)])
                                /   dx;
                        double d = ( Phi[gIndex(z,  y+1,x,  ny,nx,1)]
                                -   Phixyz)
                                /   dx;
                        double e = ( Phixyz
                                -   Phi[gIndex(z-1,y,  x,  ny,nx,1)])
                                /   dx;
                        double f = ( Phi[gIndex(z+1,y,  x,  ny,nx,1)]
                                -   Phixyz)
                                /   dx;

                        double S = Phi0xyz / sqrt(Phi0xyz*Phi0xyz + dx*dx);

                        double G = 0.0;

                        if (Phi0xyz < 0) {
                            G = std::sqrt( std::max(std::pow(std::max(a, 0.),2),std::pow(std::min(b, 0.),2))
                                    + std::max(std::pow(std::max(c, 0.),2),std::pow(std::min(d, 0.),2))
                                    /*+ max(pow(max(e, 0.),2),pow(min(f, 0.),2))*/
                                    ) - 1.;
                        } else if (Phi0xyz > 0) {
                            G = std::sqrt( std::max(std::pow(std::min(a, 0.),2),std::pow(std::max(b, 0.),2))
                                    + std::max(std::pow(std::min(c, 0.),2),std::pow(std::max(d, 0.),2))
                                    /*+ max(pow(min(e, 0.),2),pow(max(f, 0.),2))*/
                                    ) - 1.;
                        }

                        V[gIndex(z,y,x,ny,nx,0)] = S * G;
                        //if (z == nz/2 and y == 9 and x == 7) std::cout << Phi0xyz << ", " << Phixyz << ", " << S << ", " << G << ", " << -S*G << std::endl;
                    }
                }
            }

            double Vmax = 0;
            for (int i = 0; i < nx*ny*nz; ++i) {
                Vmax = std::max(Vmax, fabs(V[i]));
            }

            double dt = 0.5*dx/Vmax;

            #pragma omp parallel for collapse(3)
            for (int z = 0; z < nz; ++z) {
                for (int y = 0; y < ny; ++y) {
                    for (int x = 0; x < nx; ++x) {
                        Phi[gIndex(z,y,x,ny,nx,1)] = Phi[gIndex(z,y,x,ny,nx,1)] + dt * V[gIndex(z,y,x,ny,nx,0)];
                    }
                }
            }

            free(V);
        }

        free(Phi0);
    }

private:

    double dx;
    double numericEps;
    double maxCurvature;

    inline double delta(double x, double eps) {
        return eps / (PI*(x*x + eps*eps));
    };

    inline double heaviside(double x, double eps) {
        return 0.5 * (1.0 + 2.0/PI * atan( x / eps)); // fejl 0.5 * (1.0 + 1.0/PI * arctan( x / eps));
    };

    double areaInside(const double * const image, const double * const Phi, int nz, int ny, int nx, double eps) {

        double pixelWeight = 0.0;
        double totalWeight = 0.0;

        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    pixelWeight += image[gIndex(z,y,x,ny,nx,0)]
                                * heaviside(Phi[gIndex(z,y,x,ny,nx,1)], eps);
                    totalWeight += heaviside(Phi[gIndex(z,y,x,ny,nx,1)], eps);
                }
            }
        }

        return pixelWeight / totalWeight;
    }

    double areaOutside(const double * const image, const double * const Phi,
                      int nz, int ny, int nx, double eps) {

        double pixelWeight = 0.0;
        double totalWeight = 0.0;

        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    pixelWeight += image[gIndex(z,y,x,ny,nx,0)]
                                * (1 - heaviside(Phi[gIndex(z,y,x,ny,nx,1)], eps));
                    totalWeight += 1 - heaviside(Phi[gIndex(z,y,x,ny,nx,1)], eps);
                }
            }
        }

        return pixelWeight / totalWeight;
    }

    void applyBoundaryCondition(double * const Phi, int nz, int ny, int nx) {

        #pragma omp parallel for collapse(3)
        for (int z = -1; z < nz+1; ++z) {
            for (int y = -1; y < ny+1; ++y) {
                for (int x = -1; x < nx+1; ++x) {

                    int cx = std::min( std::max(x,0), nx-1 );
                    int cy = std::min( std::max(y,0), ny-1 );
                    int cz = std::min( std::max(z,0), nz-1 );

                    Phi[gIndex(z,y,x,ny,nx,1)] = Phi[gIndex(cz,cy,cx,ny,nx,1)];
                }
            }
        }
    }
};

#endif // CHANVESE3D_CU_H