CXX=g++
CFLAGS=-O3

SXOutputHRefSeq: SXOutputHRefSeq.o OutputHRefSeq.o
	${CXX} ${CFLAGS} -o SXOutputHRefSeq SXOutputHRefSeq.o OutputHRefSeq.o

SXOutputHRefSeq.o: SXOutputHRefSeq.cpp OutputHRefSeq.h
	${CXX} ${CFLAGS} -c -o SXOutputHRefSeq.o SXOutputHRefSeq.cpp

OutputHRefSeq.o: OutputHRefSeq.cpp OutputHRefSeq.h
	${CXX} ${CFLAGS} -c -o OutputHRefSeq.o OutputHRefSeq.cpp
clean:
	rm -f *.o

