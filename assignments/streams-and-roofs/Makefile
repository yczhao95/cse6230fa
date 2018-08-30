
CC = icc
COPTFLAGS = -O3 -xHost -qopt-report=5
OMPFLAGS = -qopenmp
CFLAGS = -g -Wall -std=c99 -fPIC
CPPFLAGS = -I$(CUDAROOT)/include
CULIBS = -L$(CUDAROOT)/lib64 -Wl,-rpath,$(CUDAROOT)/lib64 -lcudart
LIBS = $(CULIBS)
NVCC = nvcc -ccbin=icpc
NVCCFLAGS = -lineinfo -Xcompiler '-fPIC'
CUOPTFLAGS = -O3
RM = rm -f
LINKER = icpc $(OMPFLAGS)

all: cloud

%.o: %.c
	$(CC) $(CFLAGS) $(COPTFLAGS) $(CPPFLAGS) $(OMPFLAGS) -c -o $@ $<

%.link.o: %.o
	$(NVCC) -dlink -o $@ $< -lcudadevrt

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(CUOPTFLAGS) -dc -o $@ $<

verlet.o: verlet.c verlet.h cloud.h

cloud.o: cloud.c verlet.h cloud.h

cloud: verlet.o cloud.o
	$(LINKER) -o $@ $^ -Wl,-rpath,.

clean:
	$(RM) *.o cloud

.PHONY: clean
