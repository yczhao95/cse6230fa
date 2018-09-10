
MAKERULES=intel
include Makefile.$(MAKERULES)
CPPFLAGS = -I../../utils -I../../utils/tictoc
DEFINES =
RM = rm -f
LIBS = -lm
LINKER = $(CC)

all: cloud

%.o: %.c
	$(CC) -std=c99 $(CFLAGS) $(COPTFLAGS) $(CPPFLAGS) $(DEFINES) $(OMPFLAGS) -c -o $@ $<

%.link.o: %.o
	$(NVCC) -dlink -o $@ $< -lcudadevrt

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(CUOPTFLAGS) -dc -o $@ $<

verlet.o: verlet.c verlet.h cloud.h

cloud.o: cloud.c verlet.h cloud.h

cloud: verlet.o verlet_stream_and_noise.o verlet_accelerate.o cloud.o
	$(LINKER) -o $@ $^ $(LIBS)

PERF =
STREAM_N = 1000000

NP=64
DT=0.01
K=1.0
D=1.0
NT=1000

runcloud:
	$(RM) cloud cloud.o verlet.o verlet_stream_and_noise.o verlet_accelerate.o
	$(MAKE) verlet.o
	$(MAKE) cloud
	$(PERF) ./cloud $(NP) $(NT) $(DT) $(K) $(D)

clean:
	$(RM) *.o cloud

.PHONY: clean runcloud
