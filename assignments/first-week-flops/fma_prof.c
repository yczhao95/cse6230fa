
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "fma_loop.h"

int
main (int argc, char **argv)
{
  int N, T;
  int num_threads;
  float *a = NULL, b, c;

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
  if (!N) {
    a = NULL;
  }
  else {
    a = (float *) malloc (N * sizeof (float));
    if (!a) {
      printf ("Failed to allocate a\n");
      return 1;
    }
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
      a[i] = (float) i;
    }
  }

#if defined(_OPENMP)
  num_threads = omp_get_num_threads();
#else
  num_threads = 1;
#endif

#pragma omp parallel
  {
    int my_thread;
    int my_start, my_end;
    int my_N;

#if defined(_OPENMP)
    my_thread = omp_get_thread_num();
#else
    my_thread = 0;
#endif

    my_start = ((size_t) my_thread * (size_t) N) / (size_t) num_threads;
    my_end   = ((size_t) (my_thread + 1) * (size_t) N) / (size_t) num_threads;
    my_N     = my_end - my_start;
    /* execute the loop */
    fma_loop (my_N, T, &a[my_start], b, c);
  }

  printf ("%zu flops executed\n", (size_t) N * (size_t) T * 2);

  /* clean up */
  free (a);
  return 0;
}
