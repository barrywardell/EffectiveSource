effsource: main.c phis.c
	CC --std=c99 -O3 -lm -o effsource main.c phis.c

.PHONY : clean
clean :
	rm effsource