CFLAGS  = -g -O0 -Wall -std=gnu99 -I../ -I/opt/local/include
LDFLAGS = -L/opt/local/lib
LIBS = -lm -lgsl -lgslcblas

DEPS    = ../effsource.h
OBJECTS = effsource.o decompose.o

all : kerr-equatorial kerr-circular kerr-circular-m-mode kerr-equatorial-m-mode

kerr-equatorial : kerr-equatorial.o $(OBJECTS) $(DEPS)
	gcc $(CFLAGS) -o kerr-equatorial kerr-equatorial.o $(OBJECTS) $(LDFLAGS) $(LIBS)

kerr-circular : kerr-circular.o $(OBJECTS) $(DEPS)
	gcc $(CFLAGS) -o kerr-circular kerr-circular.o $(OBJECTS) $(LDFLAGS) $(LIBS)

kerr-circular-m-mode : kerr-circular-m-mode.o $(OBJECTS) $(DEPS)
	gcc $(CFLAGS) -o kerr-circular-2D-dphi2 kerr-circular-2D-dphi2.o $(OBJECTS) $(LDFLAGS) $(LIBS)

kerr-equatorial-m-mode : kerr-equatorial-m-mode.o $(OBJECTS) $(DEPS)
	gcc $(CFLAGS) -o kerr-equatorial-m-mode kerr-equatorial-m-mode $(OBJECTS) $(LDFLAGS) $(LIBS)

kerr-equatorial.o : ../kerr-equatorial.c
	gcc $(CFLAGS) -c -o kerr-equatorial.o ../kerr-equatorial.c

kerr-circular.o : ../kerr-circular.c
	gcc $(CFLAGS) -c -o kerr-circular.o ../kerr-circular.c

kerr-equatorial-m-mode.o : ../kerr-equatorial-m-mode.c
	gcc $(CFLAGS) -c -o kerr-equatorial-m-mode.o ../kerr-equatorial-m-mode.c

kerr-circular-m-mode.o : ../kerr-circular-m-mode.c
	gcc $(CFLAGS) -c -o kerr-circular-m-mode.o ../kerr-circular-m-mode.c

.PHONY : clean
clean :
	-rm -rf kerr-equatorial kerr-circular kerr-equatorial-m-mode kerr-circular-m-mode
	-rm -rf *.o *.dSYM
