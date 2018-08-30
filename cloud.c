
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
compute_hamiltonian (int Np, real_t *H,
                     const real_t *x, const real_t *y, const real_t *z,
                     const real_t *u, const real_t *v, const real_t *w)
{
  for (int i = 0; i < Np; i++) {
    H[i] = 0.5 * (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])        /* kinetic energy */
           - 1. / sqrt (x[i]*x[i] + y[i]*y[i] + z[i]*z[i]); /* potential energy */
  }
}

static int
write_step (int Np, int k, const char *basename, const real_t *x, const real_t *y, const real_t *z, const real_t *H)
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
main (int argc, char **argv)
{
  int Np, Nt, err;
  int Nint, s;
  real_t dt;
  real_t *x, *y, *z, *u, *v, *w;
  real_t *Hin, *Hout;
  const char *basename = NULL;

  if (argc != 4 && argc != 6) {
    printf ("Usage: %s NUM_POINTS NUM_STEPS DT [OUTPUT_FREQ OUTPUT_BASENAME]\n", argv[0]);
    return 1;
  }
  Np = atoi (argv[1]);
  Nt = atoi (argv[2]);
  dt = (real_t) atof (argv[3]);
  if (argc == 6) {
    Nint = atoi (argv[4]);
    basename = argv[5];
  }
  else {
    Nint = Nt;
  }

  printf ("%s, NUM_POINTS=%d, NUM_STEPS=%d, DT=%g\n", argv[0], Np, Nt, dt);

  err = safeMALLOC (Np * sizeof (real_t), &x);CHK(err);
  err = safeMALLOC (Np * sizeof (real_t), &y);CHK(err);
  err = safeMALLOC (Np * sizeof (real_t), &z);CHK(err);
  err = safeMALLOC (Np * sizeof (real_t), &u);CHK(err);
  err = safeMALLOC (Np * sizeof (real_t), &v);CHK(err);
  err = safeMALLOC (Np * sizeof (real_t), &w);CHK(err);
  err = safeMALLOC (Np * sizeof (real_t), &Hin);CHK(err);
  err = safeMALLOC (Np * sizeof (real_t), &Hout);CHK(err);

  for (int i = 0; i < Np; i++) { /* make some orbits that are visually appealing */
    real_t e, p, theta, inc, Omega, omicr;
    real_t r, x_hat, y_hat, u_hat, v_hat;
    real_t sininc, cosinc, sinOmg, cosOmg, sinomc, cosomc;
    real_t T[3][2];

    /* create a random orbit (without too much eccentricity, and not too close
     * to the center, and at a random phase in its period) */
    e = 0.5 * ((real_t) rand() / (real_t) RAND_MAX);
    p = 0.5 + ((real_t) rand() / (real_t) RAND_MAX);
    theta = 2. * PI * ((real_t) rand() / (real_t) RAND_MAX);

    r = p / (1. + e * cos (theta));
    x_hat = r * cos (theta);
    y_hat = r * sin (theta);
    u_hat = - sqrt (1. / p) * sin (theta);
    v_hat =   sqrt (1. / p) * (e + cos (theta));

    /* randomly orient the orient relative to the reference coordinates
     * https://en.wikipedia.org/wiki/Orbital_elements#Euler_angle_transformations */
    inc   = 2. * PI * ((real_t) rand() / (real_t) RAND_MAX);
    Omega = 2. * PI * ((real_t) rand() / (real_t) RAND_MAX);
    omicr = 2. * PI * ((real_t) rand() / (real_t) RAND_MAX);

    sininc = sin (inc);
    cosinc = cos (inc);
    sinOmg = sin (Omega);
    cosOmg = cos (Omega);
    sinomc = sin (omicr);
    cosomc = cos (omicr);

    T[0][0] =   cosOmg * cosomc - sinOmg * sinomc * cosinc;
    T[0][1] =   sinOmg * cosomc + cosOmg * sinomc * cosinc;
    T[1][0] = - cosOmg * sinomc - sinOmg * cosomc * cosinc;
    T[1][1] = - sinOmg * sinomc + cosOmg * cosomc * cosinc;
    T[2][0] =   cosOmg * sininc;
    T[2][1] = - sinOmg * sininc;


    x[i] = T[0][0] * x_hat + T[0][1] * y_hat;
    y[i] = T[1][0] * x_hat + T[1][1] * y_hat;
    z[i] = T[2][0] * x_hat + T[2][1] * y_hat;
    u[i] = T[0][0] * u_hat + T[0][1] * v_hat;
    v[i] = T[1][0] * u_hat + T[1][1] * v_hat;
    w[i] = T[2][0] * u_hat + T[2][1] * v_hat;
    Hin[i] = Hout[i] = 0.5 * (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])        /* kinetic energy */
                       - 1. / sqrt (x[i]*x[i] + y[i]*y[i] + z[i]*z[i]); /* potential energy */
  }
  compute_hamiltonian (Np, Hin, x, y, z, u, v, w);
  memcpy (Hout, Hin, Np * sizeof (real_t));

  for (s = 0; s < Nt; s += Nint) {
    if (basename) {err = write_step (Np, s / Nint, basename, x, y, z, Hout);CHK(err);}
    verlet_step (Np, Nint, dt, x, y, z, u, v, w);
    compute_hamiltonian (Np, Hout, x, y, z, u, v, w);
  }
  if (basename) {
    err = write_step (Np, s / Nint, basename, x, y, z, Hout);CHK(err);
  }

  printf ("Hamiltonian relative error:\n");
  for (int i = 0; i < Np; i++) {
    printf ("[%d]: %g\n", i, fabs (Hout[i] - Hin[i]) / fabs (Hin[i]));
  }

  free (Hout);
  free (Hin);
  free (w);
  free (v);
  free (u);
  free (z);
  free (y);
  free (x);
  return 0;
}
