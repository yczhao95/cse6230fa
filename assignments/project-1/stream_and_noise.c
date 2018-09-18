
#include "cloud.h"
#include "verlet.h"

void
stream_and_noise (int Np, double dt_stream, double dt_noise,
                  cse6230rand_t *rand,
                  double *restrict X[3], const double *restrict U[3])
{
  size_t tag;

  tag = cse6230rand_get_tag (rand);
  #pragma omp for schedule(static)
  for (int i = 0; i < Np; i+= 4) {
    double rval[3][4];

    for (int d = 0; d < 3; d++) {
      cse6230rand_normal_hash (rand, tag, i, d, 0, &rval[d][0]);
    }

    if (i + 4 <= Np) {
      for (int d = 0; d < 3; d++) {
        for (int j = 0; j < 4; j++) {
          X[d][i + j] += dt_stream * U[d][i + j] + dt_noise * rval[d][j];
        }
      }
    }
    else {
      for (int d = 0; d < 3; d++) {
        for (int j = 0; j < Np - i; j++) {
          X[d][i + j] += dt_stream * U[d][i + j] + dt_noise * rval[d][j];
        }
      }
    }
  }
}
