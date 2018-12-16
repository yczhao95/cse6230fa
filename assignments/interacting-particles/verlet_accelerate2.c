
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
			for (int j = i + 1; j < Np; j++) {
				double du[3];
				force (dt, X[0][i], X[1][i], X[2][i], X[0][j], X[1][j], X[2][j], du);
				for (int d = 0; d < 3; d++) {
					U_private[d][i] += du[d];
					U_private[d][j] -= du[d];
				}
			}
		}
    #pragma omp critical
		for (int i = 0; i < Np; i++) {
			for (int d = 0; d < 3; d++) {
				U[d][i] += U_private[d][i];
			}
		}
		}

}
