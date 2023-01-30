all: brain

brain: brain.c CGP-Library-V2.4/src/cgp.c CGP-Library-V2.4/src/cgp.h tests/cartpole.c tests/cartpole.h
		gcc -o brain brain.c CGP-Library-V2.4/src/cgp.c tests/cartpole.c -O3 -lm

clean: 
		-rm -f brain