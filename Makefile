CFLAGS  = -g -O0 -Wall -std=gnu99 -I/opt/local/include
LDFLAGS = -L/opt/local/lib
LIBS = -lm -lgsl -lgslcblas

DEPS    = effsource.h
OBJECTS = main.o

all : kerr-equatorial kerr-circular kerr-circular-2D

kerr-equatorial : kerr-equatorial.o $(OBJECTS) $(DEPS)
	gcc -o kerr-equatorial kerr-equatorial.o $(OBJECTS) $(LDFLAGS) $(LIBS)

kerr-circular : kerr-circular.o $(OBJECTS) $(DEPS)
	gcc -o kerr-circular kerr-circular.o $(OBJECTS) $(LDFLAGS) $(LIBS)

kerr-circular-2D : kerr-circular-2D.o $(OBJECTS) $(DEPS)
	gcc -o kerr-circular-2D kerr-circular-2D.o $(OBJECTS) $(LDFLAGS) $(LIBS)

kerr-circular-2D-dphi2 : kerr-circular-2D-dphi2.o $(OBJECTS) $(DEPS)
	gcc -o kerr-circular-2D-dphi2 kerr-circular-2D-dphi2.o $(OBJECTS) $(LDFLAGS) $(LIBS)

kerr-circular-2D-num : kerr-circular-2D.o main2.o decompose.o $(DEPS)
	gcc -o kerr-circular-2D-num kerr-circular-2D.o main2.o decompose.o $(LDFLAGS) $(LIBS)

kerr-circular-2D-num-dphi2 : kerr-circular-2D-dphi2.o main2.o decompose.o $(DEPS)
	gcc -o kerr-circular-2D-num-dphi2 kerr-circular-2D-dphi2.o main2.o decompose.o $(LDFLAGS) $(LIBS)

.PHONY : clean
clean :
	-rm -rf kerr-equatorial kerr-equatorial.dSYM kerr-circular kerr-circular.dSYM kerr-circular-2D kerr-circular-2D.dSYM
	-rm -rf $(OBJECTS) kerr-equatorial.o kerr-circular.o kerr-circular-2D.o
