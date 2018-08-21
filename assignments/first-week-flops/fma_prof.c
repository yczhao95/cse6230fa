
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "fma_omp.h"
#include "fma_cuda.h"

#define CHK(err) do {if (err) {fprintf(stderr, "[%s, %d] Top level error\n", __FILE__, __LINE__); return err;}} while (0)

int
main (int argc, char **argv)
{
  int Nh, Nd, T;
  float *ah = NULL;
  float b, c;
  float **ad = NULL;
  int numDevices = 0;
  int err;

  /* input processing */
  if (argc != 6) {
    printf ("Usage: %s HOST_ARRAY_SIZE DEV_ARRAY_SIZE LOOP_COUNT b c\n", argv[0]);
    return 1;
  }
  Nh = atoi (argv[1]);
  if (Nh < 0) {
    printf ("HOST_ARRAY_SIZE negative\n");
    return 1;
  }
  Nd = atoi (argv[2]);
  if (Nd < 0) {
    printf ("DEV_ARRAY_SIZE negative\n");
    return 1;
  }
  T = atoi (argv[3]);
  if (T < 0) {
    printf ("LOOP_COUNT negative\n");
    return 1;
  }
  b = atof (argv[4]);
  c = atof (argv[5]);

  /* initialize the array */
  err = fma_dev_initialize (Nd, T, &numDevices, &ad); CHK (err);
  err = fma_host_initialize (Nh, T, &ah); CHK (err);

  err = fma_dev_start (Nd, T, numDevices, ad, b, c); CHK (err);
  err = fma_host_start (Nh, T, ah, b, c); CHK (err);
  err = fma_dev_end (Nd, T, numDevices, ad, b, c); CHK (err);
  err = fma_host_end (Nh, T, ah, b, c); CHK (err);

  printf ("\n[%s]: %zu flops executed\n\n", argv[0], (size_t) (Nh + Nd * numDevices) * (size_t) T * 2);

  /* clean up */
  err = fma_host_free (Nh, T, &ah);  CHK (err);
  err = fma_dev_free (Nd, T, &numDevices, &ad);  CHK (err);
  return 0;
}
