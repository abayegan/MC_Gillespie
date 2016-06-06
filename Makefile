#Makefile
mc: mc.o misc.o
	gcc -o mc mc.o misc.o -L. -lRNA -lm

mc.o: mc.c misc.h misc.c
	gcc -g -c mc.c

misc.o : misc.c misc.h
	gcc -g -c misc.c

clean :
	rm -f *.o mc
