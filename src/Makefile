
CC=g++

all: tn_frnzy

tn_frnzy: dio.o matlabfunctions.o platinum.o star.o synthesis.o test.o wavread.o
	$(CC) -g -o $@ $^ -L. -lm -lfftw3

.cpp.o:
	$(CC) -I. -Wall -g -I/usr/include -c $<

clean:
	rm -f *.o

