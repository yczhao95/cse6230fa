
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "cse6230rand.h"
#include "cloud_util.h"
#include "jacobi.h"

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

int
main (int argc, char **argv)
{
  int err, size, rank;
  int p[3];  /* processes per direction */
  int t[3];  /* threads per direction */
  int q[3];  /* my rank coordinates */
  int b[3];  /* total boxes per direction */
  int bl[3]; /* local (owned) boxes per direction */
  int bh[3]; /* local + halo boxes per direction */
  int *bounds;
  int nthreads = 1;
  int nbh;   /* total number of local + halo boxes */
  int nbxy;  /* total number of local + halo z rows */
  int np;    /* total number of particles */
  int npl;   /* local number of particles */
  int nt;    /* total number of Jacobi sweeps */
  double *u, *ucopy, *f; /* arrays where solution and rhs are stored */
  double ***U, ***Ucopy, ***F; /* pointers that allow for triple index notation */
  MPI_Datatype sendtype[3][2]; /* subarrays to send, one for each face in each direction */
  MPI_Datatype recvtype[3][2]; /* subarrays to recv, one for each face in each direction */
  int          neigh[3][2];    /* neighbor for halo exchange */
  MPI_Comm comm, cartcomm;

  /* The first thing in every MPI program */
  err = MPI_Init_thread (&argc, &argv, MPI_THREAD_FUNNELED, NULL);
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

#if defined(_OPENMP)
  nthreads = omp_get_max_threads();
#endif

  /* Get the arguments */
  if (argc != 12) {
    if (!rank) { /* only the first process should print */
      printf ("Usage: %s MPI_X MPI_Y MPI_Z OMP_X OMP_Y OMP_Z B_X B_Y B_Z N_P N_T\n\n", argv[0]);
      printf ("MPI_X: processes in x direction\n");
      printf ("MPI_Y: processes in y direction\n");
      printf ("MPI_Z: processes in z direction\n");
      printf ("OMP_X: threads in x direction\n");
      printf ("OMP_Y: threads in y direction\n");
      printf ("OMP_Z: threads in z direction\n");
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

  /* create a cartesian communicator */
  {
    int periods[3] = {1, 1, 1};
    err = MPI_Cart_create (comm, 3, p, periods, 0, &cartcomm); MPI_CHK(err);
    err = MPI_Cart_coords (cartcomm, rank, 3, q); MPI_CHK(err);
  }

  /* get the number of local boxes in each direction */
  for (int i = 0; i < 3; i++) {
    int box_start;
    int box_end;

    b[i] = atoi (argv[i + 7]);

    box_start = (q[i] * b[i]) / p[i];
    box_end   = ((q[i] + 1) * b[i]) / p[i];
    bl[i] = box_end - box_start;
    bh[i] = bl[i] + 2;
    /* make the z halo 64 bytes wide and round to a multiple of 64 bytes to
     * put z rows on cache line boundaries */
    if (i == 2) {
      bh[i] += 6;
      bh[i] += 4 - (bh[i] % 4);
    }
  }

  /* get the number of threads in each direction */
  for (int i = 0; i < 3; i++) {
    t[i] = atoi (argv[i + 4]);
  }

  if (t[0] * t[1] * t[2] != nthreads) {
    if (!rank) {
      printf ("Bad thread grid: %d x %d x %d != %d\n", t[0], t[1], t[2], nthreads);
    }
    return -1;
  }

  err = safeMALLOC (16 * nthreads * sizeof (int), &bounds); CHK(err);

  compute_thread_bounds (bl, t, bounds);

  /* get the total number of particles */
  np = atoi (argv[10]);
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
  nt = atoi (argv[11]);

  if (!rank) {
    printf ("[%d x %d x %d] processes, [%d x %d x %d] threads per process, [%d x %d x %d] boxes, %d particles, %d smoother steps\n",
            p[0], p[1], p[2], t[0], t[1], t[2], b[0], b[1], b[2], np, nt);
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
      /* again the, z halo buffer is wider */
      startsend[2] += 3;
      startrecv[2] += 3;
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
  err = safeMALLOC (nbh * sizeof (double), &ucopy); CHK(err);

  /* allocate pointers that we can use to index u and f using local
   * indices, U[i][j][k], for -1 <= i <= bl[0] and so on */
  err = safeMALLOC (bh[0] * sizeof (double **), &U); CHK(err);
  err = safeMALLOC (bh[0] * sizeof (double **), &Ucopy); CHK(err);
  err = safeMALLOC (bh[0] * sizeof (double **), &F); CHK(err);
  for (int i = 0; i < bh[0]; i++) {
    err = safeMALLOC (bh[1] * sizeof (double *), &U[i]); CHK(err);
    err = safeMALLOC (bh[1] * sizeof (double *), &Ucopy[i]); CHK(err);
    err = safeMALLOC (bh[1] * sizeof (double *), &F[i]); CHK(err);
  }
  #pragma omp parallel for collapse(2) schedule(static)
  for (int i = 0; i < bh[0]; i++) {
    for (int j = 0; j < bh[1]; j++) {
      /* set U[i][j] / F[i][j] to point at the first non ghost boxes
       * in the 3D array */
      U[i][j] = &u[(i*bh[1] + j)*bh[2] + 4];
      Ucopy[i][j] = &u[(i*bh[1] + j)*bh[2] + 4];
      F[i][j] = &f[(i*bh[1] + j)*bh[2] + 4];
    }
  }
  for (int i = 0; i < bh[0]; i++) {
    (U[i])++;
    (Ucopy[i])++;
    (F[i])++;
  }
  U++;
  Ucopy++;
  F++;

  #pragma omp parallel
  {
    int tid = 0;

#if defined(_OPENMP)
    tid = omp_get_thread_num();
#endif
    for (int i = bounds[16 * tid + 0]; i < bounds[16 * tid + 1]; i++) {
      for (int j = bounds[16 * tid + 2]; j < bounds[16 * tid + 3]; j++) {
        for (int k = bounds[16 * tid + 4]; k < bounds[16 * tid + 5]; k++) {
          F[i][j][k] = 0.;
          U[i][j][k] = 0.;
        }
      }
    }
  }

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
    double srtime = 0.;

    for (int t = 0; t < nt; t++) {
      MPI_Request req[12];

      /* update halos */
      srtime -= MPI_Wtime();
      if (size > 1) {
        for (int d = 0; d < 3; d++) {
          for (int j = 0; j < 2; j++) {
            err = MPI_Isend (u, 1, sendtype[d][j], neigh[d][j], d*2 + j,
                cartcomm, &req[2 * (d*2 + j)]); MPI_CHK(err);
            err = MPI_Irecv (u, 1, recvtype[d][j^1], neigh[d][j^1], d*2 + j,
                cartcomm, &req[2 * (d*2 + j) + 1]); MPI_CHK(err);
          }
        }
        err = MPI_Waitall(12, req, MPI_STATUSES_IGNORE); MPI_CHK(err);
      }
      else {
        /* copy periodic boundaries with threads */
        #pragma omp parallel
        {
          int tid = 0;
          int bx[2], by[2], bz[2];

#if defined(_OPENMP)
          tid = omp_get_thread_num();
#endif

          bx[0] = bounds[16*tid + 0];
          bx[1] = bounds[16*tid + 1];
          by[0] = bounds[16*tid + 2];
          by[1] = bounds[16*tid + 3];
          bz[0] = bounds[16*tid + 4];
          bz[1] = bounds[16*tid + 5];
          if (bx[0] == 0) {
            for (int j = by[0]; j < by[1]; j++) {
              for (int k = bz[0]; k < bz[1]; k++) {
                U[-1][j][k] = U[bl[0]-1][j][k];
              }
            }
          }
          if (bx[1] == bl[0]) {
            for (int j = by[0]; j < by[1]; j++) {
              for (int k = bz[0]; k < bz[1]; k++) {
                U[bl[0]][j][k] = U[0][j][k];
              }
            }
          }
          if (by[0] == 0) {
            for (int i = bx[0]; i < bx[1]; i++) {
              for (int k = bz[0]; k < bz[1]; k++) {
                U[i][-1][k] = U[i][bl[1]-1][k];
              }
            }
          }
          if (by[1] == bl[1]) {
            for (int i = bx[0]; i < bx[1]; i++) {
              for (int k = bz[0]; k < bz[1]; k++) {
                U[i][bl[1]][k] = U[i][0][k];
              }
            }
          }
          if (bz[0] == 0) {
            for (int i = bx[0]; i < bx[1]; i++) {
              for (int j = by[0]; j < by[1]; j++) {
                U[i][j][-1] = U[i][j][bl[2]-1];
              }
            }
          }
          if (bz[1] == bl[1]) {
            for (int i = bx[0]; i < bx[1]; i++) {
              for (int j = by[0]; j < by[1]; j++) {
                U[i][j][bl[2]] = U[i][j][0];
              }
            }
          }
        }
      }
      srtime += MPI_Wtime();

      jacobi_sweep (bl, (const int *restrict) bounds, (const double *restrict *restrict *restrict ) F, (double *restrict *restrict *restrict) U, (double *restrict *restrict *restrict) Ucopy);
    }
    time += MPI_Wtime();

    if (!rank) {
      printf ("Sweep time: %g seconds (%g seconds communication)\n", time, srtime);
      printf ("Rate: %g lattice updates per second\n", ((double) b[0] * (double) b[1] * (double) b[2] * (double) nt) / time);
      printf ("Rate: %g halo exchanges per MPI process per second\n", (double) nt / (double) size / srtime);
    }
  }

  /* clean up */
  --U;
  --Ucopy;
  --F;
  for (int i = 0; i < bh[0]; i++) {
    free (--(U[i]));
    free (--(Ucopy[i]));
    free (--(F[i]));
  }
  free (U);
  free (Ucopy);
  free (F);
  free (u);
  free (f);
  free (bounds);

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
