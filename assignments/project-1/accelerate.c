
#include "cloud_util.h"
#include "accelerate.h"
#include "steric.h"
#include "interactions.h"

struct _accel_t
{
  int Np;
  double L;
  double k;
  double r;
  int use_ix;
  IX ix;
};

int
AccelCreate(int Np, double L, double k, double r, int use_ix, Accel *accel)
{
	int err;
  Accel a;

  err = safeMALLOC(sizeof(*a), &a);CHK(err);
  a->Np = Np;
  a->L  = L;
  a->k  = k;
  a->r  = r;
  a->use_ix = use_ix;
  if (use_ix) {
    int boxdim = L / (r * 2.); /* how could we choose boxdim ? */
		double ptl_vol = (4./3.) * 3.1415926 * r*r*r;
		double domain_vol = L*L*L;
		//maxNx will always be Np since it's possible for all particles to be existing in one box
		//an estimation would be
    int maxNx =  Np * Np * 8 * ptl_vol / domain_vol / 2;/* how should we estimate the maximum number of interactions? */
		printf("boxdim%d, maxNx%d",boxdim, maxNx);
    err = IXCreate(L, boxdim, maxNx, &(a->ix));CHK(err);
	}
  else {
    a->ix = NULL;
	}
  *accel = a;
  return 0;
}

int
AccelDestroy(Accel *accel)
{
  int err;

  if ((*accel)->ix) {
    err = IXDestroy(&((*accel)->ix));
		CHK(err);
	}
  free (*accel);
  *accel = NULL;
  return 0;
}

static void
accelerate_ix (Accel accel, Vector X, Vector U)
{
	int err;
  IX ix = accel->ix;
  int Np = X->Np;
  int Npairs;
  ix_pair *pairs;
  double L = accel->L;
  double k = accel->k;
  double r = accel->r;
  #pragma omp parallel for
  for (int i = 0; i < Np; i++) {
    for (int j = 0; j < 3; j++) {
      IDX(U,j,i) = 0.;
    }
  }

  IXGetPairs (ix, X, 2.*r, &Npairs, &pairs);	
 /* 
  #pragma omp parallel
  {
	Vector U_private;
	VectorCreate(Np, &U_private);
  #pragma omp for
  for (int i = 0; i < Np; i++) {
    for (int j = 0; j < 3; j++) {
      IDX(U_private,j,i) = 0.;
    }
  }*/	
	
  //#pragma omp parallel for
  for (int p = 0; p < Npairs; p++) {
    int i = pairs[p].p[0];
    int j = pairs[p].p[1];
    double du[3];
    force (k, r, L, IDX(X,0,i), IDX(X,1,i), IDX(X,2,i), IDX(X,0,j), IDX(X,1,j), IDX(X,2,j), du);
    for (int d = 0; d < 3; d++) {
      IDX(U,d,i) += du[d];
			IDX(U,d,j) -= du[d];
			//printf("i:%d,j:%di\n",i,j);
			//printf("value:%f, %f\n", IDX(U_private,d,i), IDX(U_private,d,j));
    }
  }
/* 
	#pragma omp critical
	for (int p = 0; p < Npairs; p++) {
    int i = pairs[p].p[0];
    int j = pairs[p].p[1];
		for (int d = 0; d < 3; d++) {
			printf("i:%d,j:%d, Np:%d\n",i,j, Np);
			printf("value:%f, %f\n", IDX(U_private,d,i), IDX(U_private,d,j));
			double pi = IDX(U_private,d,i);
			double pj = IDX(U_private,d,j);
			IDX(U,d,i) += pi;
			IDX(U,d,j) += pj;
		}
	}*/
  //VectorDestroy(&U_private);
	//}
  IXRestorePairs (ix, X, 2.*r, &Npairs, &pairs);
}

static void
accelerate_direct (Accel accel, Vector X, Vector U)
{
  int Np = accel->Np;
  double L = accel->L;
  double k = accel->k;
  double r = accel->r;

  #pragma omp parallel for schedule(static)
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

void
accelerate (Accel accel, Vector X, Vector U)
{
  if (accel->use_ix) {
    accelerate_ix (accel, X, U);
  }
  else {
    accelerate_direct (accel, X, U);
  }
}


