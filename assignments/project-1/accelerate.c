
#include "verlet.h"
#include "cloud.h"

void
accelerate (int Np, double k, double r, double L, const double *restrict X[3], double *restrict U[3])
{
  for (int i = 0; i < Np; i++) {
    double u[3] = {0.};

    for (int j = 0; j < Np; j++) {
      if (j != i) {
        double du[3];

        force (k, r, L, X[0][i], X[1][i], X[2][i], X[0][j], X[1][j], X[2][j], du);

        for (int d = 0; d < 3; d++) {
          u[d] += du[d];
        }
      }
    }
    for (int d = 0; d < 3; d++) {
      /* Instead of adding to the velocity,
       * For this project the computed interactions give the complete velocity */
      U[d][i] = u[d];
    }
  }
}

