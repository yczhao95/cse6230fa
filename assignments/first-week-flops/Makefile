
all: fma_prof fma_prof_opt

CC = icc
COPTFLAGS = -O3 -xHost -qopt-report=5
OMPFLAGS = -qopenmp
CFLAGS = -g -Wall -std=c99
CPPFLAGS =
LIBS = -l
RM = rm -f

%.o: %.c
	$(CC) $(CFLAGS) $(COPTFLAGS) $(CPPFLAGS) $(OMPFLAGS) -c -o $@ $<

fma_prof: fma_prof.o fma_loop.o
	$(CC) $(OMPFLAGS) -o $@ $^

fma_prof_opt: fma_prof.o fma_loop_opt.o
	$(CC) $(OMPFLAGS) -o $@ $^

N = $(shell echo $$(( $(PBS_NP) * 256 )))
T = 67108864
b = 0.5
c = 3.0

PERF = perf stat -v
OMPENV = OMP_PROC_BIND=spread OMP_NUM_THREADS=$(PBS_NP)

run_fma_prof: fma_prof fma_prof_opt
	$(OMPENV) $(PERF) ./fma_prof $(N) $(T) $(b) $(c)
	$(OMPENV) $(PERF) ./fma_prof_opt $(N) $(T) $(b) $(c)

clean:
	$(RM) *.o *.optrpt fma_prof fma_prof_opt

.PHONY: clean run_fma_prof
