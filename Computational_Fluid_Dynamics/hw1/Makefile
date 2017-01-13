CC = g++
CFLAGS = -O2 -fopenmp -std=c++11

derivative: main.cpp derivative.h
	$(CC) $(CFLAGS) -o hw1 main.cpp

clean:
	rm -f hw1 *~ \#* *.dat
