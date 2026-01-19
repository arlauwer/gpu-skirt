#include <cuda_runtime.h>
#include "test_gpu.cuh"

int get_gpu_device_count()
{
    int count = 0;
    // Standard CUDA runtime call to see how many GPUs are visible
    cudaError_t err = cudaGetDeviceCount(&count);

    if (err != cudaSuccess)
    {
        return -1;  // Indicates a driver or library error
    }
    return count;
}
