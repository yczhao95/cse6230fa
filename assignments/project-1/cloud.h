#if !defined(CLOUD_H)
#define      CLOUD_H
#include <math.h>
#include <cse6230rand.h>
#include "steric.h"

void initialize_variables (int Np, double k, cse6230rand_t *rand, double *X0[3], double *X[3], double *U[3]);
double compute_hamiltonian (int Np, double k, const double *X[3], const double *U[3]);

#endif
