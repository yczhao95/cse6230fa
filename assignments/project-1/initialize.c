
#include "cloud.h"

void
initialize_positions (int Np, double k, double L, cse6230rand_t *rand, double *X0[3], double *X[3])
{
  size_t init_tag;

  init_tag = cse6230rand_get_tag (rand);
  #pragma omp parallel for schedule(static)
  for (int i = 0; i < Np; i+=4) { /* for every particle */
    double xval[3][4];

    for (int d = 0; d < 3; d++) { /* get four random doubles for each variable */
      cse6230rand_hash (rand, init_tag, i,     d, 0, &xval[d][0]);
    }

    if (i + 4 <= Np) {
      for (int j = 0; j < 4; j++) {
        for (int d = 0; d < 3; d++) { /* scale uniform [0,1) variables to [-L/2, L/2) */
          X0[d][i + j] = X[d][i + j] = L * xval[d][j] - L / 2.;
        }
      }
    }
    else {
      for (int j = 0; j < Np - i; j++) {
        for (int d = 0; d < 3; d++) {
          X0[d][i + j] = X[d][i + j] = L * xval[d][j] - L / 2.;
        }
      }
    }
  }
}

