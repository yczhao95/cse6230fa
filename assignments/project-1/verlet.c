#include <math.h>
#include "cloud.h"
#include "verlet.h"

void
verlet_step (int Np, int Nt, double dt, double k, double d,
             double r, double L,
             cse6230rand_t *rand,
             double *restrict X[3],
             double *restrict U[3])
{
  int t;
  double dt_noise = sqrt(2. * d * dt);

  for (t = 0; t < Nt; t++) {
    accelerate (Np, k, r, L, (const double *restrict *) X, U);
    stream_and_noise (Np, dt, dt_noise, rand, X, (const double *restrict *) U);
  }
}

