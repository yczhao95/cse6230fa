
CC = icc
CFLAGS = -g -Wall -O3 -std=c99
CPPFLAGS =
RM = rm -f
LIBS = -lm

all: cloud

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<


verlet.o: verlet.c verlet.h cloud.h

cloud.o: cloud.c verlet.h cloud.h

cloud: verlet.o cloud.o
	$(CC) -o $@ $^ $(LIBS)

clean:
	$(RM) *.o cloud

.PHONY: clean
