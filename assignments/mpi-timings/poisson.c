
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "cse6230rand.h"
#include "cloud_util.h"

/* decode an MPI error code */
int
decode_MPI_error (int err, const char *file, const char *func, int line)
{
  char errbuf[BUFSIZ];
  int err2;
  int strlen = BUFSIZ - 1;

  err2 = MPI_Error_string (err, errbuf, &strlen);
  if (err2 == MPI_SUCCESS) {
    printf ("(%s, %s, %d), MPI error: \"%s\"\n", file, func, line, errbuf);
  }
  else {
    printf ("(%s, %s, %d), Unknown MPI error\n", file, func, line);
  }

  return err2;
}

/* macro to decode errors with line information when they occur */
#define MPI_CHK(err)                                             \
do {                                                             \
  if (err) {                                                     \
    (void) decode_MPI_error (err, __FILE__, __func__, __LINE__); \
    return err;                                                  \
  }                                                              \
} while (0)

void
jacobi_sweep (const int bl[], const double ***F, double ***U)
{
  double **Ux[3];

  for (int d = 0; d < 3; d++) {
    Ux[d] = U[d - 1];
  }
  for (int i = 0; i < bl[0]; i++) {
    double *restrict Uxy[3][3];
    const double **Fx = F[i];

    for (int d = 0; d < 3; d++) {
      for (int e = 0; e < 3; e++) {
        Uxy[d][e] = Ux[d][e - 1];
      }
    }
    for (int j = 0; j < bl[1]; j++) {
      const double *restrict Fxy = Fx[j];
      for (int k = 0; k < bl[2]; k++) {
        Uxy[1][1][k] = (1./6.) * (Fxy[k] + Uxy[0][1][k] +
                                           Uxy[2][1][k] +
                                           Uxy[1][0][k] +
                                           Uxy[1][2][k] +
                                           Uxy[1][1][k-1] +
                                           Uxy[1][1][k+1]);
      }
      for (int d = 0; d < 3; d++) {
        Uxy[d][0] = Uxy[d][1];
        Uxy[d][1] = Uxy[d][2];
        Uxy[d][2] = Ux [d][j+2];
      }
    }
    Ux[0] = Ux[1];
    Ux[1] = Ux[2];
    Ux[2] = U [i+2];
  }
}

int
main (int argc, char **argv)
{
  int err, size, rank;
  int p[3];  /* processes per direction */
  int q[3];  /* my rank coordinates */
  int b[3];  /* total boxes per direction */
  int bl[3]; /* local (owned) boxes per direction */
  int bh[3]; /* local + halo boxes per direction */
  int nbh;   /* total number of local + halo boxes */
  int nbxy;  /* total number of local + halo z rows */
  int np;    /* total number of particles */
  int npl;   /* local number of particles */
  int nt;    /* total number of Jacobi sweeps */
  int lo_bound[3]; /* lower bound on the local domain */
  double hi_bound[3]; /* higher bound on the local domain */
  double *u, *f;
  double ***U, ***F;
  MPI_Datatype sendtype[3][2]; /* subarrays to send, one for each face in each direction */
  MPI_Datatype recvtype[3][2]; /* subarrays to recv, one for each face in each direction */
  int          neigh[3][2];    /* neighbor for halo exchange */
  MPI_Comm comm, cartcomm;

  /* The first thing in every MPI program */
  err = MPI_Init (&argc, &argv);
  if (err != MPI_SUCCESS) {
    printf ("MPI initialization failed!\n");
    return err;
  }
  /* the global communicator for all processes */
  comm = MPI_COMM_WORLD;
  /* get the total number of processes */
  err = MPI_Comm_size (comm, &size); MPI_CHK(err);
  /* get my rank */
  err = MPI_Comm_rank (comm, &rank); MPI_CHK(err);

  /* Get the arguments */
  if (argc != 9) {
    if (!rank) { /* only the first process should print */
      printf ("Usage: %s P_X P_Y P_Z B_X B_Y B_Z N_P N_T\n\n", argv[0]);
      printf ("P_X: processes in x direction\n");
      printf ("P_Y: processes in y direction\n");
      printf ("P_Z: processes in z direction\n");
      printf ("B_X: boxes in x direction\n");
      printf ("B_Y: boxes in y direction\n");
      printf ("B_Z: boxes in z direction\n");
      printf ("N_P: total number of particles\n");
      printf ("N_T: total number of Jacobi smoother sweeps\n");
    }
    return -1;
  }

  /* get the number of processes in each direction */
  for (int i = 0; i < 3; i++) {
    p[i] = atoi (argv[i + 1]);
  }

  if (p[0] * p[1] * p[2] != size) {
    if (!rank) {
      printf ("Bad process grid: %d x %d x %d != %d\n", p[0], p[1], p[2], size);
    }
    return -1;
  }

  /* get the number of boxes in each direction */
  for (int i = 0; i < 3; i++) {
    b[i] = atoi (argv[i + 4]);

    if (b[i] % p[i]) {
      if (!rank) {
        printf ("Boxes/processes mismatch: %d processes does not divide %d boxes\n", p[i], b[i]);
      }
      return -1;
    }
    bl[i] = b[i] / p[i];
    bh[i] = bl[i] + 2;
  }

  /* get the total number of particles */
  np = atoi (argv[7]);
  { /* get the local number of particles */
    int my_first = (rank * np) / size;
    int my_last = ((rank + 1) * np) / size;
    int np_check;

    npl = my_last - my_first;

    err = MPI_Allreduce (&npl, &np_check, 1, MPI_INT, MPI_SUM, comm); MPI_CHK(err);
    if (np_check != np) {
      if (!rank) {
        printf ("Failed to divide particles!\n");
        return -2;
      }
    }
  }
  /* get the total number of Jacobi sweeps */
  nt = atoi (argv[8]);

  if (!rank) {
    printf ("[%d x %d x %d] processes, [%d x %d x %d] boxes, %d particles, %d smoother steps\n",
            p[0], p[1], p[2], b[0], b[1], b[2], np, nt);
  }

  /* create a cartesian communicator */
  {
    int periods[3] = {1, 1, 1};
    err = MPI_Cart_create (comm, 3, p, periods, 0, &cartcomm); MPI_CHK(err);
    err = MPI_Cart_coords (cartcomm, rank, 3, q); MPI_CHK(err);
    for (int i = 0; i < 3; i++) {
      lo_bound[i] = q[i] * bl[i];
      hi_bound[i] = (q[i] + 1) * bl[i];
    }
  }

  /* create the subarray datatypes */
  for (int i = 0; i < 3; i++) {
    int sub[3];

    for (int j = 0; j < 3; j++) {
      sub[j] = bl[j];
    }
    sub[i] = 1;

    for (int k = 0; k < 2; k++) {
      int startsend[3];
      int startrecv[3];

      for (int j = 0; j < 3; j++) {
        startsend[j] = (j == i) ? ( k ? bl[j]     : 1 ) : 1;
        startrecv[j] = (j == i) ? ( k ? bl[j] + 1 : 0 ) : 1;
      }
      err = MPI_Type_create_subarray (3, bh, sub, startsend, MPI_ORDER_C, MPI_DOUBLE, &sendtype[i][k]); MPI_CHK(err);
      err = MPI_Type_commit (&sendtype[i][k]); MPI_CHK(err);
      err = MPI_Type_create_subarray (3, bh, sub, startrecv, MPI_ORDER_C, MPI_DOUBLE, &recvtype[i][k]); MPI_CHK(err);
      err = MPI_Type_commit (&recvtype[i][k]); MPI_CHK(err);
      err = MPI_Cart_shift (cartcomm, i, k ? 1 : -1, &rank, &neigh[i][k]); MPI_CHK(err);
    }
  }

  /* allocate and zero the raw data */
  nbh = bh[0] * bh[1] * bh[2];
  err = safeMALLOC (nbh * sizeof (double), &f); CHK(err);
  err = safeMALLOC (nbh * sizeof (double), &u); CHK(err);
  for (int i = 0; i < nbh; i++) {
    f[i] = 0.;
    u[i] = 0.;
  }

  /* allocate pointers that we can use to index u and f using local
   * indices, U[i][j][k], for -1 <= i <= bl[0] and so on */
  err = safeMALLOC (bh[0] * sizeof (double **), &U); CHK(err);
  err = safeMALLOC (bh[0] * sizeof (double **), &F); CHK(err);
  for (int i = 0; i < bh[0]; i++) {
    err = safeMALLOC (bh[1] * sizeof (double *), &U[i]); CHK(err);
    err = safeMALLOC (bh[1] * sizeof (double *), &F[i]); CHK(err);
    for (int j = 0; j < bh[1]; j++) {
      /* set U[i][j] / F[i][j] to point at the first non ghost boxes
       * in the 3D array */
      U[i][j] = &u[(i*bh[1] + j)*bh[2] + 1];
      F[i][j] = &f[(i*bh[1] + j)*bh[2] + 1];
    }
    (U[i])++;
    (F[i])++;
  }
  U++;
  F++;

  /* deposit random charges */
  {
    int seed = 0;
    size_t tag  = 0;
    cse6230rand_t rand;

    cse6230rand_seed (seed, &rand);
    for (int i = 0; i < np; i += 4) {
      double rval[3][4];
      int b[3];

      /* create random coordinates */
      for (int d = 0; d < 3; d++) {
        cse6230rand_hash (&rand, tag, (size_t) rank, (size_t) i, (size_t) d, &rval[d][0]);
      }

      /* scale the local coordinates so that they belong to this process's
       * subdomain, and figure out which local box they belong to */
      if (i + 4 <= np) {
        for (int k = 0; k < 4; k++) {
          for (int j = 0; j < 3; j++) {
            b[j] = rval[j][k] * bl[j];
          }
          F[b[0]][b[1]][b[2]] += 1.;
        }
      }
      else {
        for (int k = 0; k < np - i; k++) {
          for (int j = 0; j < 3; j++) {
            b[j] = rval[j][k] * bl[j];
          }
          F[b[0]][b[1]][b[2]] += 1.;
        }
      }
    }
  }

  /* jacobi sweeps */
  {
    double time = -MPI_Wtime();

    for (int t = 0; t < nt; t++) {
      /* update halos */
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
          err = MPI_Sendrecv (u, 1, sendtype[i][j],   neigh[i][j],   i * 2 + j,
                              u, 1, recvtype[i][j^1], neigh[i][j^1], i * 2 + j,
                              cartcomm, MPI_STATUS_IGNORE); MPI_CHK(err);
        }
      }
      jacobi_sweep (bl, (const double ***) F, U);
    }
    time += MPI_Wtime();

    if (!rank) {
      printf ("Sweep time: %g seconds\n", time);
      printf ("Rate: %g lattice updates per second\n", (b[0] * b[1] * b[2] * nt) / time);
    }
  }

  /* clean up */
  --U;
  --F;
  for (int i = 0; i < bh[0]; i++) {
    free (--(U[i]));
    free (--(F[i]));
  }
  free (U);
  free (F);
  free (u);
  free (f);

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      err = MPI_Type_free (&sendtype[i][j]); MPI_CHK(err);
      err = MPI_Type_free (&recvtype[i][j]); MPI_CHK(err);
    }
  }

  err = MPI_Comm_free (&cartcomm); MPI_CHK(err);
  err = MPI_Finalize();
  return err;
}
