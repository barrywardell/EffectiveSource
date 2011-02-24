effsource: main.c phis.c phisr9.c src.c srcr9.c
	CC --std=c99 -O3 -lm -o effsource main.c phis.c phisr9.c src.c srcr9.c

.PHONY : clean
clean :
	rm effsource