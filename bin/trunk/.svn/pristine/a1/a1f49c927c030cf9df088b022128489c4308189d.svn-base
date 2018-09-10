CXX=g++

CFLAGS=-Wall -O3 -Werror
LFLAGS=-lpthread

SXReadSNV:	main.o CSXReadSNV.o SXDBSNPgChecker.o SXDBIndelgChecker.o SXAnnotate.o SXDBIndelgChecker_OdbSnpfmt.o
	${CXX} ${CFLAGS} ${LFLAGS} -o SXReadSNV_CG main.o CSXReadSNV.o SXDBSNPgChecker.o SXDBIndelgChecker.o SXAnnotate.o SXDBIndelgChecker_OdbSnpfmt.o

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

SXDBIndelgChecker_OdbSnpfmt.o: SXDBIndelgChecker_OdbSnpfmt.cpp SXDBIndelgChecker_OdbSnpfmt.h
	${CXX} ${CFLAGS} -c -o SXDBIndelgChecker_OdbSnpfmt.o SXDBIndelgChecker_OdbSnpfmt.cpp
clean:
	rm -f *.o
