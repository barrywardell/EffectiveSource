CFLAGS  = -g -O3 -Wall -std=c99
LDFLAGS = -lm

DEPS    = effsource.h
OBJECTS = main.o

all : kerr-equatorial kerr-circular

kerr-equatorial : kerr-equatorial.o $(OBJECTS) $(DEPS)
	gcc -o kerr-equatorial kerr-equatorial.o $(OBJECTS)

kerr-circular : kerr-circular.o $(OBJECTS) $(DEPS)
	gcc -o kerr-circular kerr-circular.o $(OBJECTS)

.PHONY : clean
clean :
	-rm -rf kerr-equatorial kerr-equatorial.dSYM kerr-circular kerr-circular.dSYM
	-rm -rf $(OBJECTS) kerr-equatorial.o kerr-circular.o
