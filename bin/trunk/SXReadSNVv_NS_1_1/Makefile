CXX=g++

CFLAGS=-Wall -O3 -Werror
LFLAGS=-lpthread

SX ReadSNV:	main.o CSXReadSNV.o SXDBSNPgChecker.o SXDBIndelgChecker.o SXAnnotate.o
	${CXX} ${CFLAGS} ${LFLAGS} -o SXReadSNV main.o CSXReadSNV.o SXDBSNPgChecker.o SXDBIndelgChecker.o SXAnnotate.o

main.o:	main.cpp CSXReadSNV.h
	${CXX} ${CFLAGS} -c -o main.o main.cpp

CSXReadSNV.o: CSXReadSNV.cpp CSXReadSNV.h typedef.h SXDBSNPgChecker.h SXDBIndelgChecker.h SXAnnotate.h
	${CXX} ${CFLAGS} -c -o CSXReadSNV.o CSXReadSNV.cpp

SXDBSNPgChecker.o: SXDBSNPgChecker.cpp SXDBSNPgChecker.h
	${CXX} ${CFLAGS} -c -o SXDBSNPgChecker.o SXDBSNPgChecker.cpp

SXDBIndelgChecker.o: SXDBIndelgChecker.cpp SXDBIndelgChecker.h
	${CXX} ${CFLAGS} -c -o SXDBIndelgChecker.o SXDBIndelgChecker.cpp

SXAnnotate.o: SXAnnotate.cpp SXAnnotate.h
	${CXX} ${CFLAGS} -c -o SXAnnotate.o SXAnnotate.cpp
clean:
	rm -f *.o
