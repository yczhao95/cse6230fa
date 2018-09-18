#if !defined(CLOUD_H)
#define      CLOUD_H
#include <math.h>
#include <cse6230rand.h>
#include "steric.h"

void initialize_positions (int Np, double k, double L, cse6230rand_t *rand, double *X0[3], double *X[3]);

#endif
