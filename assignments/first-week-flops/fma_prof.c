
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
  int N, T;
  float *ah = NULL;
  float b, c;
  float **ad = NULL;
  int numDevices = 0;
  int err;

  /* input processing */
  if (argc != 5) {
    printf ("Usage: %s ARRAY_SIZE LOOP_COUNT b c\n", argv[0]);
    return 1;
  }
  N = atoi (argv[1]);
  if (N < 0) {
    printf ("ARRAY_SIZE negative\n");
    return 1;
  }
  T = atoi (argv[2]);
  if (T < 0) {
    printf ("LOOP_COUNT negative\n");
    return 1;
  }
  b = atof (argv[3]);
  c = atof (argv[4]);

  /* initialize the array */
  err = fma_dev_initialize (N, T, &numDevices, &ad); CHK (err);
  err = fma_host_initialize (N, T, &ah); CHK (err);

  err = fma_dev_start (N, T, numDevices, ad, b, c); CHK (err);
  err = fma_host_start (N, T, ah, b, c); CHK (err);
  err = fma_dev_end (N, T, numDevices, ad, b, c); CHK (err);
  err = fma_host_end (N, T, ah, b, c); CHK (err);


  printf ("\n[%s]: %zu flops executed\n\n", argv[0], (size_t) N * (size_t) T * 2);

  /* clean up */
  err = fma_host_free (N, T, &ah);  CHK (err);
  err = fma_dev_free (N, T, &numDevices, &ad);  CHK (err);
  return 0;
}
