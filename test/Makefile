CFLAGS  = -g -O0 -Wall -std=gnu99 -I../ -I/opt/local/include
LDFLAGS = -L/opt/local/lib
LIBS = -lm -lgsl -lgslcblas

DEPS    = ../effsource.h
OBJECTS = effsource.o decompose.o

all : kerr-equatorial kerr-circular

kerr-equatorial : kerr-equatorial.o kerr-equatorial-coeffs.o kerr-equatorial-dtcoeffs.o kerr-equatorial-dttcoeffs.o $(OBJECTS) $(DEPS)
	gcc $(CFLAGS) -o kerr-equatorial kerr-equatorial-coeffs.o kerr-equatorial-dtcoeffs.o kerr-equatorial-dttcoeffs.o kerr-equatorial.o $(OBJECTS) $(LDFLAGS) $(LIBS)

kerr-circular : kerr-circular.o $(OBJECTS) $(DEPS)
	gcc $(CFLAGS) -o kerr-circular kerr-circular.o $(OBJECTS) $(LDFLAGS) $(LIBS)

kerr-equatorial.o : ../kerr-equatorial.c
	gcc $(CFLAGS) -c -o kerr-equatorial.o ../kerr-equatorial.c

kerr-equatorial-coeffs.o : ../kerr-equatorial-coeffs.c
	gcc $(CFLAGS) -c -o kerr-equatorial-coeffs.o ../kerr-equatorial-coeffs.c

kerr-equatorial-dtcoeffs.o : ../kerr-equatorial-dtcoeffs.c
	gcc $(CFLAGS) -c -o kerr-equatorial-dtcoeffs.o ../kerr-equatorial-dtcoeffs.c

kerr-equatorial-dttcoeffs.o : ../kerr-equatorial-dttcoeffs.c
	gcc $(CFLAGS) -c -o kerr-equatorial-dttcoeffs.o ../kerr-equatorial-dttcoeffs.c

kerr-circular.o : ../kerr-circular.c
	gcc $(CFLAGS) -c -o kerr-circular.o ../kerr-circular.c

.PHONY : clean
clean :
	-rm -rf kerr-equatorial kerr-circular
	-rm -rf *.o *.dSYM
