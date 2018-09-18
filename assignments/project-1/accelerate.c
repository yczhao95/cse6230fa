
#include "cloud_util.h"
#include "accelerate.h"
#include "steric.h"

struct _accel_t
{
  int Np;
  double L;
  double k;
  double r;
};

int
AccelCreate(int Np, double L, double k, double r, Accel *accel)
{
  int err;
  Accel a;

  err = safeMALLOC(sizeof(*a), &a);CHK(err);
  a->Np = Np;
  a->L  = L;
  a->k  = k;
  a->r  = r;
  *accel = a;
  return 0;
}

int
AccelDestroy(Accel *accel)
{
  free (*accel);
  *accel = NULL;
  return 0;
}

void
accelerate (Accel accel, Vector X, Vector U)
{
  int Np = accel->Np;
  double L = accel->L;
  double k = accel->k;
  double r = accel->r;

  #pragma omp for schedule(static)
  for (int i = 0; i < Np; i++) {
    double u[3] = {0.};

    for (int j = 0; j < Np; j++) {
      if (j != i) {
        double du[3];

        force (k, r, L, IDX(X,0,i), IDX(X,1,i), IDX(X,2,i), IDX(X,0,j), IDX(X,1,j), IDX(X,2,j), du);

        for (int d = 0; d < 3; d++) {
          u[d] += du[d];
        }
      }
    }
    for (int d = 0; d < 3; d++) {
      /* Instead of adding to the velocity,
       * For this project the computed interactions give the complete velocity */
      IDX(U,d,i) = u[d];
    }
  }
}

