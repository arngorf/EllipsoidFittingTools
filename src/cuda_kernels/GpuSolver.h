#ifndef GPUSOLVER_H
#define GPUSOLVER_H

/*
 * Example code
 */

class GpuInterface {
public:
    int n[20];
    int y;
    int asize;

    GpuInterface();
    int calculateSum();
    void setY(int);
};

#endif // GPUSOLVER_H