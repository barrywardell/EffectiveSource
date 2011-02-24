effsource: main.c effsource.c
	CC -Wall --std=c99 -O3 -lm -o effsource main.c effsource.c

.PHONY : clean
clean :
	rm effsource