
all: fma_prof fma_prof_opt

CC = icc
COPTFLAGS = -O3 -xHost -qopt-report=5
OMPFLAGS = -qopenmp
CFLAGS = -g -Wall -std=c99 -fPIC
CPPFLAGS = -I$(CUDAROOT)/include
CULIBS = -L$(CUDAROOT)/lib64 -lcudart
NVCC = nvcc -ccbin=icpc -cudart=static
NVCCFLAGS = -G -Xcompiler '-fPIC'
NVCCOPTFLAGS = 
RM = rm -f
LINKER = icpc $(OMPFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $(COPTFLAGS) $(CPPFLAGS) $(OMPFLAGS) -c -o $@ $<

%.link.o: %.o
	$(NVCC) -dlink -o $@ $< -lcudadevrt

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(NVCCOPTFLAGS) -dc -o $@ $<

fma_cuda_link.o: fma_cuda.o fma_loop_dev.o
	$(NVCC) $(NVCCFLAGS) -dlink  $^ -o $@

libfma_cuda.so: fma_cuda_link.o fma_cuda.o fma_loop_dev.o
	$(LINKER) -shared -Wl,-soname,libfma_cuda.so -o $@ $^ $(CULIBS)

libfma_cuda_opt.so: fma_cuda.o fma_loop_dev_opt.o
	$(NVCC) -lib $^ -o $@

fma_prof: fma_prof.o fma_omp.o fma_loop_host.o libfma_cuda.so
	$(LINKER) -o $@ fma_prof.o fma_omp.o fma_loop_host.o -L. -ldl -Wl,-rpath,. -lfma_cuda $(CULIBS) 

fma_prof_opt: fma_prof.o fma_omp.o fma_loop_host_opt.o libfma_cuda_opt.so
	$(LINKER) -o $@ $^ $(LIBS)

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
	$(RM) *.o *.optrpt *.so fma_prof fma_prof_opt

.PHONY: clean run_fma_prof
