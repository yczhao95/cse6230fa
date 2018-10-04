
#include <stdio.h>
#if defined(_OPENMP)
#include <omp.h>
#endif

void
compute_thread_bounds (const int *restrict bl, const int *restrict t, int *restrict bounds)
{
  #pragma omp parallel
  {
    int tid = 0;
    int bstart[3];
    int bend[3];
    int tidx, tidy, tidyx, tidz;

#if defined(_OPENMP)
    tid = omp_get_thread_num();
#endif

    tidz = tid % t[2];
    tidyx = tid / t[2];
    tidy = tidyx % t[1];
    tidx = tidyx / t[1];

    bounds[16*tid + 0] = (tidx * bl[0]) / t[0];
    bounds[16*tid + 1] = ((tidx + 1) * bl[0]) / t[0];
    bounds[16*tid + 2] = (tidy * bl[1]) / t[1];
    bounds[16*tid + 3] = ((tidy + 1) * bl[1]) / t[1];
    bounds[16*tid + 4] = (tidz * bl[2]) / t[2];
    bounds[16*tid + 5] = ((tidz + 1) * bl[2]) / t[2];

#if 0
    #pragma omp critical
    {
      printf ("%d (%d, %d, %d): [%d, %d] x [%d, %d] x, [%d, %d]\n",
              tid, tidx, tidy, tidz,
              bounds[16*tid+0],
              bounds[16*tid+1],
              bounds[16*tid+2],
              bounds[16*tid+3],
              bounds[16*tid+4],
              bounds[16*tid+5]);
      fflush(stdout);
    }
#endif
  }
}

void
jacobi_sweep (const int *restrict bl, const int *restrict bounds,
              const double *restrict *restrict *restrict F,
              double *restrict *restrict *restrict U,
              double *restrict *restrict *restrict Ucopy)
{
  #pragma omp parallel
  {
    int tid = 0;
    int bstart[3];
    int bend[3];

#if defined(_OPENMP)
    tid = omp_get_thread_num();
#endif

    bstart[0] = bounds[16 * tid + 0];
    bend[0]   = bounds[16 * tid + 1];
    bstart[1] = bounds[16 * tid + 2];
    bend[1]   = bounds[16 * tid + 3];
    bstart[2] = bounds[16 * tid + 4];
    bend[2]   = bounds[16 * tid + 5];
    for (int i = bstart[0]; i < bend[0]; i++) {
      for (int j = bstart[1]; j < bend[1]; j++) {
        double *restrict Uxy[3][3];
        double *restrict uc = Ucopy[i][j];
        const double *restrict f = F[i][j];

        for (int d = 0; d < 3; d++) {
          for (int e = 0; e < 3; e++) {
            Uxy[d][e] = U[i + d - 1][j + e - 1];
          }
        }
        for (int k = bstart[2]; k < bend[2]; k++) {
          uc[k] = Uxy[1][1][k] + (1./6.) * (f[k] + Uxy[0][1][k] +
                                                   Uxy[2][1][k] +
                                                   Uxy[1][0][k] +
                                                   Uxy[1][2][k] +
                                                   Uxy[1][1][k-1] +
                                                   Uxy[1][1][k+1]);
        }
      }
    }
    #pragma omp barrier
    for (int i = bstart[0]; i < bend[0]; i++) {
      for (int j = bstart[1]; j < bend[1]; j++) {
        const double *restrict uc = Ucopy[i][j];
        double *restrict u = U[i][j];

        for (int k = bstart[2]; k < bend[2]; k++) {
          u[k] = uc[k];
        }
      }
    }
  }
}
