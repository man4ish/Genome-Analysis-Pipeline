CXX=g++
CFLAGS=-g

SXGetMMRefSeq: SXGetMMRefSeq.o GetMMRefSeq.o
	${CXX} ${CFLAGS} -o SXGetMMRefSeq SXGetMMRefSeq.o GetMMRefSeq.o

SXGetMMRefSeq.o: SXGetMMRefSeq.cpp GetMMRefSeq.h
	${CXX} ${CFLAGS} -c -o SXGetMMRefSeq.o SXGetMMRefSeq.cpp

GetMMRefSeq.o: GetMMRefSeq.cpp GetMMRefSeq.h
	${CXX} ${CFLAGS} -c -o GetMMRefSeq.o GetMMRefSeq.cpp
clean:
	rm -f *.o

