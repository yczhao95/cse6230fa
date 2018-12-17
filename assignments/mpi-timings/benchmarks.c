
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define MPI_CHK(err) if (err != MPI_SUCCESS) return err

#define MPI_LOG(rank,...) if (!rank) {printf(__VA_ARGS__);fflush(stdout);}
#define max(a,b) ((a) > (b) ? (a) : (b))
int splitCommunicator(MPI_Comm comm, int firstCommSize, MPI_Comm *subComm_p)
{
  /* TODO: split the communicator `comm` into one communicator for ranks
   * [0, ..., firstCommSize - 1] and one for [firstCommSize, ..., size - 1],
   * where `size` is the size of `comm`.
   *
   * Look at the way subcommunicators are created using a coloring in
   * `VecGetGlobalSize_tree_subcomm()` in the file `dmv_global_size.c` in the
   * `notes/mpi/dmv` example.
   */
	int err, rank, color;
	err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);
	color = rank >= firstCommSize;
	err = MPI_Comm_split(comm,color,firstCommSize,subComm_p); MPI_CHK(err);
  return 0;
}

int destroyCommunicator(MPI_Comm *subComm_p)
{
  /* TODO: destroy the subcommunicator created in `splitCommunicator` */
	int err;
	err = MPI_Comm_free(subComm_p); MPI_CHK(err);
  return 0;
}

int startTime(MPI_Comm comm, double *tic_p)
{
  /* TODO: insert a barrier on `comm` to synchronize the processes */
  /* TODO: Record the MPI walltime in `tic_p` */
	MPI_Barrier(comm);
	*tic_p = MPI_Wtime();
  return 0;
}

int stopTime(double tic_in, double *toc_p)
{
  /* TODO: Get the elapsed MPI walltime since `tic_in`,
   * write the results in `toc_p` */
	*toc_p = MPI_Wtime() - tic_in;
  return 0;
}

int maxTime(MPI_Comm comm, double myTime, double *maxTime_p)
{
  /* TODO: take the times from all processes and compute the maximum,
   * storing the result on process 0 */
	*maxTime_p = max(*maxTime_p, myTime); 
  return 0;
}

int main(int argc, char **argv)
{
  MPI_Comm comm;
  int      err;
  int      size, rank;
  int      maxSize = 1 << 24;
  int      maxCollectiveSize;
  int      numTests = 10;
  char     *buffer;
  char     *buffer2;

  err = MPI_Init(&argc, &argv); MPI_CHK(err);

  comm = MPI_COMM_WORLD;

  err = MPI_Comm_size(comm, &size); MPI_CHK(err);
  err = MPI_Comm_rank(comm, &rank); MPI_CHK(err);

  maxCollectiveSize = maxSize / size;

  /* Print out some info about the environment */
  if (!rank) {
    char library[MPI_MAX_LIBRARY_VERSION_STRING] = {0};
    int version, subversion, liblen;
    double time, timeprec;

    err = MPI_Get_version(&version, &subversion); MPI_CHK(err);
    err = MPI_Get_library_version(library, &liblen); MPI_CHK(err);

    printf("MPI Version: %d.%d\n", version, subversion);
    printf("%s\n", library);
    printf("MPI # Procs: %d\n", size);

    time = MPI_Wtime();
    timeprec = MPI_Wtick();

    printf("MPI Wtime %g, precision %g\n", time, timeprec);
    if (MPI_WTIME_IS_GLOBAL) {
      printf("MPI Wtime is global\n");
    } else {
      printf("MPI Wtime is not global\n");
    }
  }
  for (int i = 0; i < size; i++) {
    char procname[MPI_MAX_PROCESSOR_NAME] = {0};
    if (i == rank) {
      int namelen;

      err = MPI_Get_processor_name(procname, &namelen); MPI_CHK(err);
      if (i) {
        err = MPI_Send(procname, namelen, MPI_CHAR, 0, 0, comm); MPI_CHK(err);
      }
    }
    if (!rank) {
      if (i) {
        err = MPI_Recv(procname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, i, 0, comm, MPI_STATUS_IGNORE); MPI_CHK(err);
      }
      printf("MPI proc %d host: %s\n", i, procname);
    }
  }

  buffer = (char *) calloc(maxSize, sizeof(char));
  if (!buffer) return 1;
  buffer2 = (char *) calloc(maxCollectiveSize, sizeof(char));
  if (!buffer2) return 1;

  /* === POINT TO POINT TIMINGS === */
  MPI_LOG(rank, "MPI Point-to-point test:\n"
          " # Processes  | Message Size | Total Size   | Time         | B/s\n");
  for (int numComm = 4; numComm <= size; numComm *= 4) {
    for (int numBytes = 8; numBytes <= maxSize; numBytes *= 8) {
      double        timeAvg = 0.;
      long long int totalNumBytes = numBytes * numComm;

      for (int t = 0; t < numTests; t++) {
        double tic = -1.;

        err = startTime(comm, &tic); MPI_CHK(err);
        if (rank < numComm) {
          if (rank % 2) {
            err = MPI_Send(buffer, numBytes, MPI_BYTE, rank - 1, 0, comm); MPI_CHK(err);
            err = MPI_Recv(buffer, numBytes, MPI_BYTE, rank - 1, 0, comm, MPI_STATUS_IGNORE); MPI_CHK(err);
          } else {
            err = MPI_Recv (buffer, numBytes, MPI_BYTE, rank + 1, 0, comm, MPI_STATUS_IGNORE); MPI_CHK(err);
            err = MPI_Send (buffer, numBytes, MPI_BYTE, rank + 1, 0, comm); MPI_CHK(err);
          }
        }
        err = stopTime(tic, &tic); MPI_CHK(err);
        if (t) {
          timeAvg += tic;
        }
      }
      timeAvg /= (numTests - 1);
      MPI_LOG(rank, " %12d   %12d   %12lld   %+12.5e   %+12.5e\n", numComm, numBytes, totalNumBytes, timeAvg, totalNumBytes / timeAvg);
    }
  }

  /* === ONE TO ALL: Proc 0 broadcasts messages of varying lengths to all other processes === */
  MPI_LOG(rank, "MPI One-to-all broadcast test:\n"
          " # Processes  | Message Size | Total Size   | Time         | B/s\n");
  for (int numComm = 4; numComm <= size; numComm *= 4) {
    MPI_Comm subComm = MPI_COMM_NULL;

    err = splitCommunicator(MPI_COMM_WORLD, numComm, &subComm); MPI_CHK(err);
    for (int numBytes = 8; numBytes <= maxSize; numBytes *= 8) {
      double        timeAvg = 0.;
      long long int totalNumBytes = numBytes * (numComm - 1);

      for (int t = 0; t < numTests; t++) {
        double tic = -1.;

        err = startTime(comm, &tic); MPI_CHK(err);
        if (rank < numComm) {
          /* TODO: Use the subComm communicator to broadcast from rank 0 the
           * first `numBytes` bytes of the `buffer` to the other subComm
           * processes. Store the results in the first `numBytes` bytes of the
           * `buffer` on the receiving processes. */
				//	MPI_Bcast(buffer, numComm, MPI_BYTE, 0,comm);
        }
        err = stopTime(tic, &tic); MPI_CHK(err);
        if (t) {
          timeAvg += tic;
        }
      }
      timeAvg /= (numTests - 1);
      err = maxTime(MPI_COMM_WORLD, timeAvg, &timeAvg); MPI_CHK(err);
      MPI_LOG(rank, " %12d   %12d   %12lld   %+12.5e   %+12.5e\n", numComm, numBytes, totalNumBytes, timeAvg, totalNumBytes / timeAvg);
    }
    err = destroyCommunicator(&subComm); MPI_CHK(err);
  }

  /* === ONE TO ALL: Proc 0 scatters messages of varying lengths to all other processes === */
  MPI_LOG(rank, "MPI One-to-all scatter test:\n"
          " # Processes  | Message Size | Total Size   | Time         | B/s\n");
  for (int numComm = 4; numComm <= size; numComm *= 4) {
    MPI_Comm subComm = MPI_COMM_NULL;

    err = splitCommunicator(MPI_COMM_WORLD, numComm, &subComm); MPI_CHK(err);
    for (int numBytes = 8; numBytes <= maxCollectiveSize; numBytes *= 8) {
      double        timeAvg = 0.;
      long long int totalNumBytes = numBytes * (numComm - 1);

      for (int t = 0; t < numTests; t++) {
        double tic = -1.;

        err = startTime(comm, &tic); MPI_CHK(err);
        if (rank < numComm) {
          /* TODO: Use the subComm communicator to scatter from rank 0 the first
           * (`numBytes` x `size`) bytes of `buffer` to the other processes
           * (each process receives `numBytes` worth). Store the results in the
           * first `numBytes` bytes of `buffer2` on the receiving processes.
           */
        }
        err = stopTime(tic, &tic); MPI_CHK(err);
        if (t) {
          timeAvg += tic;
        }
      }
      timeAvg /= (numTests - 1);
      err = maxTime(MPI_COMM_WORLD, timeAvg, &timeAvg); MPI_CHK(err);
      MPI_LOG(rank, " %12d   %12d   %12lld   %+12.5e   %+12.5e\n", numComm, numBytes, totalNumBytes, timeAvg, totalNumBytes / timeAvg);
    }
    err = destroyCommunicator(&subComm); MPI_CHK(err);
  }

  /* === ALL TO ONE: Proc 0 minimizes messages of varying lengths from all other processes === */
  MPI_LOG(rank, "MPI All-to-one reduce test:\n"
          " # Processes  | Message Size | Total Size   | Time         | B/s\n");
  for (int numComm = 4; numComm <= size; numComm *= 4) {
    MPI_Comm subComm = MPI_COMM_NULL;

    err = splitCommunicator(MPI_COMM_WORLD, numComm, &subComm); MPI_CHK(err);
    for (int numBytes = 8; numBytes <= maxSize; numBytes *= 8) {
      double        timeAvg = 0.;
      long long int totalNumBytes = numBytes * (numComm - 1);

      for (int t = 0; t < numTests; t++) {
        double tic = -1.;

        err = startTime(comm, &tic); MPI_CHK(err);
        if (rank < numComm) {
          /* TODO: Use the subComm communicator to compute the minimum of the
           * first `numBytes` chars of `buffer` from every process and store the
           * results in `buffer` on process 0 (HINT: look up the proper usage of
           * MPI_IN_PLACE).
           */
        }
        err = stopTime(tic, &tic); MPI_CHK(err);
        if (t) {
          timeAvg += tic;
        }
      }
      timeAvg /= (numTests - 1);
      err = maxTime(MPI_COMM_WORLD, timeAvg, &timeAvg); MPI_CHK(err);
      MPI_LOG(rank, " %12d   %12d   %12lld   %+12.5e   %+12.5e\n", numComm, numBytes, totalNumBytes, timeAvg, totalNumBytes / timeAvg);
    }
    err = destroyCommunicator(&subComm); MPI_CHK(err);
  }

  /* === ALL TO ONE: Proc 0 gathers messages of varying lengths from all other processes === */
  MPI_LOG(rank, "MPI All-to-one gather test:\n"
          " # Processes  | Message Size | Total Size   | Time         | B/s\n");
  for (int numComm = 4; numComm <= size; numComm *= 4) {
    MPI_Comm subComm = MPI_COMM_NULL;

    err = splitCommunicator(MPI_COMM_WORLD, numComm, &subComm); MPI_CHK(err);
    for (int numBytes = 8; numBytes <= maxCollectiveSize; numBytes *= 8) {
      double        timeAvg = 0.;
      long long int totalNumBytes = numBytes * (numComm - 1);

      for (int t = 0; t < numTests; t++) {
        double tic = -1.;

        err = startTime(comm, &tic); MPI_CHK(err);
        if (rank < numComm) {
          /* TODO: Use the subComm communicator to gather the first `numBytes`
           * bytes of `buffer2` from every process and store the results in
           * `buffer` on process 0.
           */
        }
        err = stopTime(tic, &tic); MPI_CHK(err);
        if (t) {
          timeAvg += tic;
        }
      }
      timeAvg /= (numTests - 1);
      err = maxTime(MPI_COMM_WORLD, timeAvg, &timeAvg); MPI_CHK(err);
      MPI_LOG(rank, " %12d   %12d   %12lld   %+12.5e   %+12.5e\n", numComm, numBytes, totalNumBytes, timeAvg, totalNumBytes / timeAvg);
    }
    err = destroyCommunicator(&subComm); MPI_CHK(err);
  }

  /* === ALL TO ALL: All processes sum messages of varying lengths === */
  MPI_LOG(rank, "MPI All-to-all reduce test:\n"
          " # Processes  | Message Size | Total Size   | Time         | B/s\n");
  for (int numComm = 4; numComm <= size; numComm *= 4) {
    MPI_Comm subComm = MPI_COMM_NULL;

    err = splitCommunicator(MPI_COMM_WORLD, numComm, &subComm); MPI_CHK(err);
    for (int numBytes = 8; numBytes <= maxSize; numBytes *= 8) {
      double        timeAvg = 0.;
      long long int totalNumBytes = numBytes * ((numComm - 1) + (numComm - 1));

      for (int t = 0; t < numTests; t++) {
        double tic = -1.;

        err = startTime(comm, &tic); MPI_CHK(err);
        if (rank < numComm) {
          /* TODO: Use the subComm communicator to sum the first `numBytes`
           * chars of `buffer` from every process and store the results in
           * `buffer` on all processes (HINT: MPI_IN_PLACE again).
           */
        }
        err = stopTime(tic, &tic); MPI_CHK(err);
        if (t) {
          timeAvg += tic;
        }
      }
      timeAvg /= (numTests - 1);
      err = maxTime(MPI_COMM_WORLD, timeAvg, &timeAvg); MPI_CHK(err);
      MPI_LOG(rank, " %12d   %12d   %12lld   %+12.5e   %+12.5e\n", numComm, numBytes, totalNumBytes, timeAvg, totalNumBytes / timeAvg);
    }
    err = destroyCommunicator(&subComm); MPI_CHK(err);
  }

  /* === ALL TO ALL: All processes gather messages of varying lengths */
  MPI_LOG(rank, "MPI All-to-all gather test:\n"
          " # Processes  | Message Size | Total Size   | Time         | B/s\n");
  for (int numComm = 4; numComm <= size; numComm *= 4) {
    MPI_Comm subComm = MPI_COMM_NULL;

    err = splitCommunicator(MPI_COMM_WORLD, numComm, &subComm); MPI_CHK(err);
    for (int numBytes = 8; numBytes <= maxCollectiveSize; numBytes *= 8) {
      double        timeAvg = 0.;
      long long int totalNumBytes = numBytes * ((numComm - 1) * numComm);

      for (int t = 0; t < numTests; t++) {
        double tic = -1.;

        err = startTime(comm, &tic); MPI_CHK(err);
        if (rank < numComm) {
          /* TODO: Use the subComm communicator to gather the first `numBytes`
           * bytes of `buffer2` from every process and store the results in
           * `buffer`.
           */
        }
        err = stopTime(tic, &tic); MPI_CHK(err);
        if (t) {
          timeAvg += tic;
        }
      }
      timeAvg /= (numTests - 1);
      err = maxTime(MPI_COMM_WORLD, timeAvg, &timeAvg); MPI_CHK(err);
      MPI_LOG(rank, " %12d   %12d   %12lld   %+12.5e   %+12.5e\n", numComm, numBytes, totalNumBytes, timeAvg, totalNumBytes / timeAvg);
    }
    err = destroyCommunicator(&subComm); MPI_CHK(err);
  }

  /* === ALL TO ALL: All processes transpose messages of varying lengths using MPI_Alltoall() */
  MPI_LOG(rank, "MPI All-to-all transpose test:\n"
          " # Processes  | Message Size | Total Size   | Time         | B/s\n");
  for (int numComm = 4; numComm <= size; numComm *= 4) {
    MPI_Comm subComm = MPI_COMM_NULL;

    err = splitCommunicator(MPI_COMM_WORLD, numComm, &subComm); MPI_CHK(err);
    for (int numBytes = 8; numBytes <= maxCollectiveSize; numBytes *= 8) {
      double        timeAvg = 0.;
      long long int totalNumBytes = numBytes * ((numComm - 1) * (numComm - 1));

      for (int t = 0; t < numTests; t++) {
        double tic = -1.;

        err = startTime(comm, &tic); MPI_CHK(err);
        if (rank < numComm) {
          /* TODO: Use the subComm communicator to transpose the first
           * `numComm` * `numBytes` bytes of `buffer` from every process and
           * store the results in `buffer`.  This is another place where
           * MPI_IN_PLACE is relevant.
           */
        }
        err = stopTime(tic, &tic); MPI_CHK(err);
        if (t) {
          timeAvg += tic;
        }
      }
      timeAvg /= (numTests - 1);
      err = maxTime(MPI_COMM_WORLD, timeAvg, &timeAvg); MPI_CHK(err);
      MPI_LOG(rank, " %12d   %12d   %12lld   %+12.5e   %+12.5e\n", numComm, numBytes, totalNumBytes, timeAvg, totalNumBytes / timeAvg);
    }
    err = destroyCommunicator(&subComm); MPI_CHK(err);
  }

  free (buffer2);
  free (buffer);
  err = MPI_Finalize();
  return err;
}
