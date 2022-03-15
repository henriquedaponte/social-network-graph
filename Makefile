CXX=g++

CXXFLAGS= -std=c++11 -Wall -g

BINARIES=PA2

PA2: PA2.o main.o
	${CXX} $^ -o $@

clean:
	/bin/rm -f ${BINARIES} *.o
