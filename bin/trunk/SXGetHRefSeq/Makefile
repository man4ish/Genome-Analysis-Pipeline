CXX=g++
CFLAGS=-O3

SXGetHRefSeq: SXGetHRefSeq.o GetHRefSeq.o
	${CXX} ${CFLAGS} -o SXGetHRefSeq SXGetHRefSeq.o GetHRefSeq.o

SXGetHRefSeq.o: SXGetHRefSeq.cpp GetHRefSeq.h
	${CXX} ${CFLAGS} -c -o SXGetHRefSeq.o SXGetHRefSeq.cpp

GetHRefSeq.o: GetHRefSeq.cpp GetHRefSeq.h
	${CXX} ${CFLAGS} -c -o GetHRefSeq.o GetHRefSeq.cpp
clean:
	rm -f *.o

