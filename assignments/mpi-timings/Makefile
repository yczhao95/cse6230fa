
MPICC    = mpicc
CFLAGS   = -g -Wall -O3 -std=c99
CPPFLAGS = 
LIBS     = 

MPIRUN    = mpirun
N         = 4
MPIRUNCMD = $(MPIRUN) -np $(N)

RM       = rm -f

%.o: %.c
	$(MPICC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

benchmarks: benchmarks.o
	$(MPICC) -o $@ $^ $(LIBS)

runbenchmarks: benchmarks
	$(MPIRUNCMD) $(abspath $<)

clean:
	$(RM) *.o benchmarks

.PHONY: clean runbenchmarks
