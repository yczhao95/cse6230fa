
#include "verlet.h"
#include "cloud.h"

void
verlet_step_accelerate (int Np, double dt, const double *restrict X[3], double *restrict U[3])
{
  #pragma omp parallel
	{
	double U_private[3][Np];
	#pragma omp for
	for (int i = 0; i < Np; i++) {
		double u1 = 0., u2 = 0., u3 = 0.;
		for (int j = 0; j < Np; j++) {
			if (j != i) {
				double du[3];
				force (dt, X[0][i], X[1][i], X[2][i], X[0][j], X[1][j], X[2][j], du);

					u1 += du[0];
					u2 += du[1];
					u3 += du[2];
			} 
		}
		U_private[0][i] += u1;
		U_private[1][i] += u2;
		U_private[2][i] += u3;
	}
  #pragma omp critical
	for (int i = 0; i < Np; i++) {
		for (int d = 0; d < 3; d++) {
			U[d][i] += U_private[d][i];
		}
	}
	}
}

