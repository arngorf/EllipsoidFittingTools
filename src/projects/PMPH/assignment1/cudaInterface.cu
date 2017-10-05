#include <iostream>
#include "cuda.h"
#include <sys/time.h>
#include <time.h>
int timeval_subtract(
    struct timeval* result, struct timeval* t2,struct timeval* t1) {
    unsigned int resolution=1000000;
    long int diff = (t2->tv_usec + resolution * t2->tv_sec) -
    (t1->tv_usec + resolution * t1->tv_sec) ;
    result->tv_sec = diff / resolution;
    result->tv_usec = diff % resolution;
    return (diff<0);
}

__global__
void functionKernel(float* d_in, float *d_out) {
    const unsigned int lid = threadIdx.x; // local id inside a block
    const unsigned int gid = blockIdx.x*blockDim.x + lid; // global id
    float val = d_in[gid]/(d_in[gid]-2.3);
    d_out[gid] = val*val*val;
}

class CudaInterface {
public:

    CudaInterface();

    void serialVsParallelTest();
private:

    int n;
    int memSize;

    float* inputArray;
    float* serialResult;
    float* parallelResult;
};

CudaInterface::CudaInterface() {

    // Test parameters
    n = 753411;

    // Store of memory
    memSize = n * sizeof(float);

    // Allocate CPU memory
    inputArray = (float*) malloc(memSize);
    serialResult = (float*) malloc(memSize);
    parallelResult = (float*) malloc(memSize);

    // Initialize input array
    for (int i = 1; i <= n; ++i) {
        inputArray[i] = (float)i;
    }
}

void CudaInterface::serialVsParallelTest() {

    // Allocate GPU memory
    float* d_in;
    float* d_out;

    cudaMalloc((void**)&d_in, memSize);
    cudaMalloc((void**)&d_out, memSize);

    // copy host memory to device
    cudaMemcpy(d_in, inputArray, memSize, cudaMemcpyHostToDevice);
    // execute the kernel

    unsigned int num_threads = n;
    unsigned int mem_size = num_threads*sizeof(float);
    unsigned int block_size = 256;
    unsigned int num_blocks = ((num_threads + (block_size - 1)) / block_size);


    // Begin timing
    unsigned long int elapsed;

    struct timeval t_start_p, t_end_p, t_diff_p;
    gettimeofday(&t_start_p, NULL);

    functionKernel<<<num_blocks, block_size>>>(d_in, d_out);

    gettimeofday(&t_end_p, NULL);
    timeval_subtract(&t_diff_p, &t_end_p, &t_start_p);
    elapsed = t_diff_p.tv_sec*1e6+t_diff_p.tv_usec;
    printf("Parallel execution took %d microseconds (%.2fms)\n",elapsed,elapsed/1000.0);

    //functionKernel<<<1, n>>>(d_in, d_out);
    // copy result from device to host
    cudaMemcpy(parallelResult, d_out, memSize, cudaMemcpyDeviceToHost);
    // print result

    cudaFree(d_in);
    cudaFree(d_out);

    struct timeval t_start_s, t_end_s, t_diff_s;
    gettimeofday(&t_start_s, NULL);

    for (int i = 0; i < n; ++i) {
        float val = inputArray[i]/(inputArray[i]-2.3);
        serialResult[i] = val*val*val;
    }

    gettimeofday(&t_end_s, NULL);
    timeval_subtract(&t_diff_s, &t_end_s, &t_start_s);
    elapsed = t_diff_s.tv_sec*1e6+t_diff_s.tv_usec;
    printf("Parallel execution took %d microseconds (%.2fms)\n",elapsed,elapsed/1000.0);

    for (int i = 0; i < 100; ++i) {
        std::cout << inputArray[i] << ", " << serialResult[i] << ", " << parallelResult[i] << std::endl;
    }

    free(inputArray);
    free(serialResult);
    free(parallelResult);

}

void cudaTest() {
    CudaInterface cudaInterface;
    cudaInterface.serialVsParallelTest();
}

int main() {

    cudaTest();

    return 0;
}