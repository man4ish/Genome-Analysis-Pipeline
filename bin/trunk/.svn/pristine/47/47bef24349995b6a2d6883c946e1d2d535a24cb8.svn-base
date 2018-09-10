CXX=g++
CFLAGS=-Wall -O3 -g0 -DNDEBUG
LFLAGS=

SSNpParser:	SSNpParser.o 
	${CXX} ${CFLAGS} ${LFLAGS} -o  SSNpParser SSNpParser.o 

SSNpParser.o:	SSNpParser.cpp SSNpParser.h
	${CXX} ${CFLAGS} -c -o SSNpParser.o SSNpParser.cpp

clean:
	rm -f *.o
