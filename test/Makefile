CFLAGS  = -g -Wall -std=gnu99 -I../ -I/usr/local/include
LDFLAGS = -L/usr/local/lib
LIBS = -lm -lgsl -lgslcblas

DEPS    = ../effsource.h

all : kerr-equatorial kerr-circular kerr-circular-comoving

kerr-equatorial : kerr-equatorial.o kerr-equatorial-coeffs.o kerr-equatorial-dtcoeffs.o kerr-equatorial-dttcoeffs.o effsource.o decompose.o $(DEPS)
	gcc $(CFLAGS) -o kerr-equatorial kerr-equatorial-coeffs.o kerr-equatorial-dtcoeffs.o kerr-equatorial-dttcoeffs.o kerr-equatorial.o effsource.o decompose.o $(LDFLAGS) $(LIBS)

kerr-circular : kerr-circular.o effsource.o decompose.o $(DEPS)
	gcc $(CFLAGS) -o kerr-circular kerr-circular.o effsource.o decompose.o $(LDFLAGS) $(LIBS)

kerr-circular-comoving : kerr-circular-comoving.o effsource.o decompose.o $(DEPS)
	gcc $(CFLAGS) -o kerr-circular-comoving kerr-circular-comoving.o effsource.o decompose.o $(LDFLAGS) $(LIBS)

effsource.o : effsource.c
	gcc $(CFLAGS) -O3 -c -o effsource.o effsource.c

decompose.o : decompose.c
	gcc $(CFLAGS) -O3 -c -o decompose.o decompose.c

kerr-equatorial.o : ../kerr-equatorial.c
	gcc $(CFLAGS) -O3 -c -o kerr-equatorial.o ../kerr-equatorial.c

kerr-equatorial-coeffs.o : ../kerr-equatorial-coeffs.c
	gcc $(CFLAGS) -O0 -c -o kerr-equatorial-coeffs.o ../kerr-equatorial-coeffs.c

kerr-equatorial-dtcoeffs.o : ../kerr-equatorial-dtcoeffs.c
	gcc $(CFLAGS) -O0 -c -o kerr-equatorial-dtcoeffs.o ../kerr-equatorial-dtcoeffs.c

kerr-equatorial-dttcoeffs.o : ../kerr-equatorial-dttcoeffs.c
	gcc $(CFLAGS) -O0 -c -o kerr-equatorial-dttcoeffs.o ../kerr-equatorial-dttcoeffs.c

kerr-circular.o : ../kerr-circular.c
	gcc $(CFLAGS) -O3 -c -o kerr-circular.o ../kerr-circular.c

kerr-circular-comoving.o : ../kerr-circular-comoving.c
	gcc $(CFLAGS) -O3 -c -o kerr-circular-comoving.o ../kerr-circular-comoving.c

.PHONY : clean
clean :
	-rm -rf kerr-equatorial kerr-circular
	-rm -rf *.o *.dSYM
