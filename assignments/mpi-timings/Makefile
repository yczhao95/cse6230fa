
OMPFLAGS = -qopenmp
MPICC    = mpicc $(OMPFLAGS)
CFLAGS   = -g -Wall -O3 -xHost -std=c99
CPPFLAGS = -I../../utils
LIBS     =

MPIRUN    = mpiexec
N         = ${PBS_NP}
MPIRUNCMD = $(MPIRUN) -bind-to core

RM       = rm -f

all: benchmarks poisson

%.o: %.c
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

benchmarks: benchmarks.o
	$(MPICC) -o $@ $^ $(LIBS)

poisson: poisson.o jacobi.o cloud_util.o
	$(MPICC) -o $@ $^ $(LIBS)

runbenchmarks: benchmarks
	$(MPIRUNCMD) -n $(N) ./benchmarks | xvfb-run python figgen.py

MPI_X=1
MPI_Y=1
MPI_Z=1
OMP_X=1
OMP_Y=1
OMP_Z=1
BX=216
BY=216
BZ=216
NP=1
NT=1000
MPI_N=1
OMP_N=1

runpoisson: poisson
	OMP_NUM_THREADS=$(OMP_N) OMP_PROC_BIND=spread $(MPIRUNCMD) -n $(MPI_N) ./poisson $(MPI_X) $(MPI_Y) $(MPI_Z) $(OMP_X) $(OMP_Y) $(OMP_Z) $(BX) $(BY) $(BZ) $(NP) $(NT)

runpoissonomp: poisson
	OMP_NUM_THREADS=$(OMP_N) OMP_PROC_BIND=true ./poisson 1 1 1 $(OMP_X) $(OMP_Y) $(OMP_Z) $(BX) $(BY) $(BZ) $(NP) $(NT)

clean:
	$(RM) *.o benchmarks poisson

.PHONY: clean runbenchmarks runpoisson
