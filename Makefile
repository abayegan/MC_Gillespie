CXX             = gcc
CXXFLAGS        = -fPIC -Iheader
LDFLAGS        = -L. -lm -lRNA
#Makefile
mc: mc.o misc.o
	${CXX} ${CXXFLAGS} -o mc mc.o misc.o ${LDFLAGS}

mc.o: mc.c misc.h misc.c
	${CXX} ${CXXFLAGS} -c mc.c

misc.o : misc.c misc.h
	${CXX} ${CXXFLAGS} -c misc.c

clean :
	rm -f *.o mc
