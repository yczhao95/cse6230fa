
extern "C" {
  #include <stdio.h>
  #include "fma_cuda.h"
}
#include "fma_dev.h"

#define CUDA_CHK(cerr) do {cudaError_t _cerr = (cerr); if ((_cerr) != cudaSuccess) {fprintf(stderr,"[%s, %d] Cuda error %s\n", __FILE__, __LINE__, cudaGetErrorString(_cerr)); return 1;}} while(0)

__global__
static void
fma_initialize (int N, float *a)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int gridSize = gridDim.x * blockDim.x;

  for (int i = tid; i < N; i+= gridSize) {
    a[i] = i;
  }
}

int
fma_dev_initialize (int N, int T, int *numDevices, float ***a)
{
  float **aa = NULL;
  cudaError_t cerr;

  cerr = cudaGetDeviceCount (numDevices); CUDA_CHK(cerr);
  printf ("numDevices %d\n", *numDevices);
  if (*numDevices) {
    aa = (float **) malloc (*numDevices * sizeof (float *));
    if (!aa) {
      fprintf (stderr, "Failed to allocate aa\n");
      return 1;
    }
    for (int i = 0; i < *numDevices; i++) {
      struct cudaDeviceProp prop;
      int block, grid;

      cerr = cudaSetDevice(i); CUDA_CHK(cerr);
      cerr = cudaMalloc (&aa[i], N * sizeof (float)); CUDA_CHK(cerr);
      cerr = cudaGetDeviceProperties (&prop, i); CUDA_CHK(cerr);
      block = prop.maxThreadsPerBlock;
      grid = (N + block - 1) / block;
      fma_initialize<<<grid, block>>>(N, aa[i]);
      cerr = cudaDeviceSynchronize(); CUDA_CHK(cerr);
    }
  }
  *a = aa;
  return 0;
}

int
fma_dev_free (int N, int T, int *numDevices, float ***a)
{
  cudaError_t cerr;

  for (int i = 0; i < *numDevices; i++) {
    cerr = cudaSetDevice(i); CUDA_CHK(cerr);
    cerr = cudaFree ((*a)[i]); CUDA_CHK(cerr);
  }
  free (*a);
  *a = NULL;
  return 0;
}

int
fma_dev_start (int N, int T, int numDevices, float **a, float b, float c)
{
  cudaError_t cerr;

  for (int i = 0; i < numDevices; i++) {
    struct cudaDeviceProp prop;
    int block, grid;

    cerr = cudaSetDevice(i); CUDA_CHK(cerr);
    cerr = cudaGetDeviceProperties (&prop, i); CUDA_CHK(cerr);
    block = prop.maxThreadsPerBlock;
    grid = (N + block - 1) / block;
    fma_loop_dev<<<grid, block>>>(N, T, a[i], b, c);
  }
  return 0;
}

int
fma_dev_end (int N, int T, int numDevices, float **a, float b, float c)
{
  cudaError_t cerr;

  for (int i = 0; i < numDevices; i++) {
    cerr = cudaSetDevice(i); CUDA_CHK(cerr);
    cerr = cudaDeviceSynchronize(); CUDA_CHK(cerr);
  }
  return 0;
}
