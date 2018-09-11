
MAKERULES=intel
include Makefile.$(MAKERULES)
CPPFLAGS = -I../../utils -I../../utils/tictoc
DEFINES =
RM = rm -f
LIBS = -lm -lrt
LINKER = $(CC)

all: cloud cloud2

%.o: %.c
	$(CC) -std=c99 $(CFLAGS) $(COPTFLAGS) $(CPPFLAGS) $(DEFINES) $(OMPFLAGS) -c -o $@ $<

%.link.o: %.o
	$(NVCC) -dlink -o $@ $< -lcudadevrt

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(CUOPTFLAGS) -dc -o $@ $<

verlet.o: verlet.c verlet.h cloud.h

verlet_stream_and_noise.o: verlet_stream_and_noise.c verlet.h cloud.h

verlet_accelerate.o: verlet_accelerate.c verlet.h cloud.h

verlet_accelerate2.o: verlet_accelerate2.c verlet.h cloud.h

cloud.o: cloud.c verlet.h cloud.h

cloud: verlet.o verlet_stream_and_noise.o verlet_accelerate.o cloud.o
	$(LINKER) -o $@ $^ $(LIBS)

cloud2: verlet.o verlet_stream_and_noise.o verlet_accelerate2.o cloud.o
	$(LINKER) -o $@ $^ $(LIBS)

PERF =
STREAM_N = 1000000

NP=64
DT=0.01
K=1.0
D=1.0
NT=1000
CHUNK=-1
VIZNAME=cloud

runcloud:
	$(MAKE) clean
	$(MAKE) cloud
	$(PERF) ./cloud $(NP) $(NT) $(DT) $(K) $(D)

runcloud2:
	$(MAKE) clean
	$(MAKE) cloud2
	$(PERF) ./cloud2 $(NP) $(NT) $(DT) $(K) $(D)

vizcloud:
	$(MAKE) clean
	$(MAKE) cloud
	./cloud $(NP) $(NT) $(DT) $(K) $(D) $(CHUNK) $(VIZNAME) | python3 cloud.py

clean:
	$(RM) *.o cloud cloud2

.PHONY: clean runcloud runcloud2 vizcloud
