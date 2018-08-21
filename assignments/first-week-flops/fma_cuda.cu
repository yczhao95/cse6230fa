
extern "C" {
  #include "fma_cuda.h"
}
#include "fma_dev.h"

#define CUDA_CHK(cerr) do {cudaError_t _cerr = (cerr); if ((_cerr) != cudaSuccess) fprintf(stderr,"[%s, %d] Cuda error %s\n", __FILE__, __LINE__, cudaGetErrorString(_cerr)); return 1;} while(0)

int
fma_dev_initialize (int N, int T, int *numDevices, float ***a)
{
  *numDevices = 0;
  *a = NULL;
  return 0;
}

int
fma_dev_free (int N, int T, int *numDevices, float ***a)
{
  *numDevices = 0;
  *a = NULL;
  return 0;
}

int
fma_dev_start (int N, int T, int numDevices, float **a, float b, float c)
{
  return 0;
}

int
fma_dev_end (int N, int T, int numDevices, float **a, float b, float c)
{
  return 0;
}
