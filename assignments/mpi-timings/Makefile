
MPICC    = mpicc
CFLAGS   = -g -Wall -O3 -std=c99
CPPFLAGS = -I../../utils
LIBS     =

MPIRUN    = mpirun
N         = 4
MPIRUNCMD = $(MPIRUN) -np $(N)

RM       = rm -f

all: benchmarks poisson

%.o: %.c
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

benchmarks: benchmarks.o
	$(MPICC) -o $@ $^ $(LIBS)

poisson: poisson.o cloud_util.o
	$(MPICC) -o $@ $^ $(LIBS)

runbenchmarks: benchmarks
	$(MPIRUNCMD) $(abspath $<)

clean:
	$(RM) *.o benchmarks poisson

.PHONY: clean runbenchmarks
