CC = g++
MPI = mpic++
CFLAGS = -O2 -fopenmp -std=c++11

sort: sort.cpp bitonic.h bitonic_MPI.h
	$(MPI) $(CFLAGS) -o sort sort.cpp

clean:
	rm -f sort *~ \#* 
