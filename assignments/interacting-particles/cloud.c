
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tictoc.h>
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
write_step (int Np, int k, const char *basename, const double *x, const double *y, const double *z, const double *H)
{
  char outname[BUFSIZ];
  FILE *fp;

  snprintf (outname, BUFSIZ-1, "%s_%d.vtk", basename, k);
  fp = fopen (outname, "w");
  if (!fp) {
    fprintf (stderr, "unable to open %s for output\n", outname);
    return 1;
  }
  fprintf (fp, "# vtk DataFile Version 2.0\n");
  fprintf (fp, "Point cloud example\n");
  fprintf (fp, "ASCII\n");
  fprintf (fp, "DATASET POLYDATA\n");
  fprintf (fp, "POINTS %d FLOAT\n", Np);
  for (int i = 0; i < Np; i++) {
    fprintf (fp, "%f %f %f\n", (float) x[i], (float) y[i], (float) z[i]);
  }
  fprintf (fp, "\nPOINT_DATA %d\n", Np);
  fprintf (fp, "SCALARS Hamiltonian float\n");
  fprintf (fp, "LOOKUP_TABLE default\n");
  for (int i = 0; i < Np; i++) {
    fprintf (fp, "%f\n", (float) H[i]);
  }
  fclose (fp);
  return 0;
}

int
process_options (int argc, char **argv, int *Np, int *Nt, int *Nint, double *dt, double *k, double *d, const char **basename)
{
  if (argc < 6 || argc > 8) {
    printf ("Usage: %s NUM_POINTS NUM_STEPS DT K D [CHUNK_SIZE OUTPUT]\n", argv[0]);
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
      *basename = argv[7];
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
  const char *basename = NULL;
  cse6230rand_t rand;
  TicTocTimer time_loop;

  err = process_options (argc, argv, &Np, &Nt, &Nint, &dt, &k, &d, &basename);CHK(err);

  printf ("[%s] NUM_POINTS=%d, NUM_STEPS=%d, CHUNK_SIZE=%d, DT=%g, K=%g, D=%g\n", argv[0], Np, Nt, Nint, dt, k, d);

  for (int d = 0; d < 3; d++) {
    err = safeMALLOC (Np * sizeof (double), &X0[d]);CHK(err);
    err = safeMALLOC (Np * sizeof (double), &X[d]);CHK(err);
    err = safeMALLOC (Np * sizeof (double), &U[d]);CHK(err);
  }

  cse6230rand_seed (seed, &rand);

  initialize_variables (Np, k, &rand, X0, X, U);

  Hin = compute_hamiltonian (Np, k, (const double **)X, (const double **)U);
  printf ("[%s] Hamiltonian, T = 0: %g\n", argv[0], Hin);

  time_loop = tic();
  for (int t = 0; t < Nt; t += Nint) {
#if 0
    if (basename) {write_step (Np, t / Nint, basename, x, y, z, Hout);}
#endif

    /* execute the loop */
    verlet_step (Np, Nint, dt, k, d, &rand, X, U);
  }
  toc(&time_loop);
  Hout = compute_hamiltonian (Np, k, (const double **)X, (const double **)U);
  printf ("[%s] Hamiltonian, T = %g: %g, Relative Error: %g\n", argv[0], Nt * dt, Hout, fabs(Hout - Hin) / Hin);
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

    printf ("[%s] Average Distance Traveled: %g\n", argv[0], avgDist);
  }
#if 0
  if (basename) {write_step (Np, Nt / Nint, basename, x, y, z, Hout);}
#endif
  for (int d = 0; d < 3; d++) {
    free (U[d]);
    free (X[d]);
    free (X0[d]);
  }
  return 0;
}
