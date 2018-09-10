
#include <tictoc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cse6230rand.h>
#include "cloud.h"
#include "verlet.h"

#define PI 3.1415926535897932384626433832795L

int safe_malloc (size_t count, void *out, const char * file, const char * fn, int line)
{
  if (!count) {
    *((void **) out) = NULL;
    return 0;
  }
  *((void **) out) = malloc (count);
  if (!(*((void **) out))) {
    fprintf (stderr, "%s, %s (%d): failed to malloc %zu bytes\n", file, fn, line, count);
    return 1;
  }
  return 0;
}

#define safeMALLOC(count,out) safe_malloc (count, out, __FILE__, __func__, __LINE__)
#define CHK(Q) if (Q) return Q

static void
initialize_variables (int Np, double k, cse6230rand_t *rand, double *X0[3], double *X[3], double *U[3])
{
  size_t init_tag;

  init_tag = cse6230rand_get_tag (rand);
  for (int i = 0; i < Np; i+=4) { /* for every particle */
    double xval[3][4];
    double uval[3][4];

    for (int d = 0; d < 3; d++) { /* get four random doubles for each variable */
      cse6230rand_hash (rand, init_tag, i,     d, 0, &xval[d][0]);
      if (k) {
        cse6230rand_hash (rand, init_tag, i, 3 + d, 0, &uval[d][0]);
      }
      else {
        for (int j = 0; j < 4; j++) {
          uval[d][j] = 0.5;
        }
      }
    }

    if (i + 4 <= Np) {
      for (int j = 0; j < 4; j++) {
        for (int d = 0; d < 3; d++) { /* scale uniform [0,1) variables to [-1, 1) */
          X0[d][i + j] = X[d][i + j] = 2. * xval[d][j] - 1.;
          U[d][i + j] = 2. * uval[d][j] - 1.;
        }
      }
    }
    else {
      for (int j = 0; j < Np - i; j++) {
        for (int d = 0; d < 3; d++) {
          X0[d][i + j] = X[d][i + j] = 2. * xval[d][j] - 1.;
          U[d][i + j] = 2. * uval[d][j] - 1.;
        }
      }
    }
  }
}

static double
compute_hamiltonian (int Np, double k,
                     const double *X[3],
                     const double *U[3])
{
  double h = 0.;

  for (int i = 0; i < Np; i++) { /* for every particle */
    h += 0.5 * (U[0][i]*U[0][i] + U[1][i]*U[1][i] + U[2][i]*U[2][i]); /* kinetic energy */
    if (k) {
      for (int j = i + 1; j < Np; j++) { /* for every other particle */
        if (j != i) {
          double hj = potential (k, X[0][i], X[1][i], X[2][i], X[0][j], X[1][j], X[2][j]);

          h += hj;
        }
      }
    }
  }
  return h;
}

static int
write_step (int Np, double T, const double *X[3])
{
  int line_size = 4;

  printf ("{ \"num_points\": %d,\n", Np);
  printf ("  \"time\": %g,\n", T);
  printf("  \"X\": [\n");
  for (int d=0; d < 3; d++) {
    printf("    [\n");
    for (int l = 0; l < Np; l+= line_size) {
      printf("     ");
      if (l + line_size < Np) {
        for (int j = 0; j < line_size; j++) {
          printf(" %+9.2e,", X[d][l + j]);
        }
      }
      else {
        for (int j = 0; j < Np - 1 - l; j++) {
          printf(" %+9.2e,", X[d][l + j]);
        }
        printf(" %+9.2e", X[d][Np - 1]);
      }
      printf("\n");
    }
    printf("    ],\n");
  }
  printf("  ]\n");
  printf ("}\n");
  return 0;
}

int
process_options (int argc, char **argv, int *Np, int *Nt, int *Nint, double *dt, double *k, double *d, int *pipe)
{
  if (argc < 6 || argc > 8) {
    printf ("Usage: %s NUM_POINTS NUM_STEPS DT K D [CHUNK_SIZE PIPE_JSON]\n", argv[0]);
    return 1;
  }
  *Np = atoi (argv[1]);
  *Nt = atoi (argv[2]);
  *dt = (double) atof (argv[3]);
  *k  = (double) atof (argv[4]);
  *d  = (double) atof (argv[5]);
  if (argc > 6) {
    *Nint = atoi (argv[6]);
    if (argc == 8) {
      *pipe = atoi (argv[7]);
    }
    else {
      *pipe = 0;
    }
  }
  else {
    *Nint = *Nt;
  }
  return 0;
}


int
main (int argc, char **argv)
{
  int Np, Nt, err;
  int Nint;
  double dt;
  double k;
  double d;
  double *X0[3];
  double *X[3], *U[3];
  double Hin, Hout;
  int seed = 6230;
  int pipe = 0;
  cse6230rand_t rand;
  TicTocTimer loop_timer;
  double loop_time;

  err = process_options (argc, argv, &Np, &Nt, &Nint, &dt, &k, &d, &pipe);CHK(err);

  for (int d = 0; d < 3; d++) {
    err = safeMALLOC (Np * sizeof (double), &X0[d]);CHK(err);
    err = safeMALLOC (Np * sizeof (double), &X[d]);CHK(err);
    err = safeMALLOC (Np * sizeof (double), &U[d]);CHK(err);
  }

  cse6230rand_seed (seed, &rand);

  initialize_variables (Np, k, &rand, X0, X, U);

  Hin = compute_hamiltonian (Np, k, (const double **)X, (const double **)U);
  if (!pipe) {
    printf ("[%s] NUM_POINTS=%d, NUM_STEPS=%d, CHUNK_SIZE=%d, DT=%g, K=%g, D=%g\n", argv[0], Np, Nt, Nint, dt, k, d);
    printf ("[%s] Hamiltonian, T = 0: %g\n", argv[0], Hin);
  }
  else {
    printf ("{ \"num_points\": %d, \"k\": %e, \"d\": %e, \"dt\": %e, \"num_steps\": %d, \"step_chunk\": %d, \"hamiltonian_0\": %e }\n",
            Np, k, d, dt, Nt, Nint, Hin);
  }


  loop_timer = tic();
  for (int t = 0; t < Nt; t += Nint) {
    if (pipe) {write_step (Np, t * dt, (const double **)X);}
    /* execute the loop */
    verlet_step (Np, Nint, dt, k, d, &rand, X, U);
  }
  loop_time = toc(&loop_timer);
  if (pipe) {write_step (Np, Nt * dt, (const double **)X);}
  Hout = compute_hamiltonian (Np, k, (const double **)X, (const double **)U);
  {
    double avgDist = 0.;

    for (int i = 0; i < Np; i++) {
      double this_dist;

      this_dist = 0.;
      for (int d = 0; d < 3; d++) {
        this_dist += (X[d][i] - X0[d][i]) * (X[d][i] - X0[d][i]);
      }
      this_dist = sqrt (this_dist);
      avgDist += this_dist;
    }
    avgDist /= Np;

    if (!pipe) {
      printf ("[%s] Simulation time: %g\n", argv[0], loop_time);
      printf ("[%s] Hamiltonian, T = %g: %g, Relative Error: %g\n", argv[0], Nt * dt, Hout, fabs(Hout - Hin) / Hin);
      printf ("[%s] Average Distance Traveled: %g\n", argv[0], avgDist);
    } else {
      printf ("{ \"avg_dist\": %e, \"hamiltonian_T\": %e, \"hamiltionian_relerr\": %e, \"sim_time\": %e }\n",
              avgDist, Hout, fabs(Hout - Hin) / Hin, loop_time);
    }
  }

  for (int d = 0; d < 3; d++) {
    free (U[d]);
    free (X[d]);
    free (X0[d]);
  }
  return 0;
}
