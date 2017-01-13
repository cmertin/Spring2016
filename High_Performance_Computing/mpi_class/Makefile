CC = g++
MPI = mpic++
CFLAGS = -O2 -fopenmp -std=c++11

ex1: ex1.cpp
	$(MPI) $(CFLAGS) -o ex1 ex1.cpp


qsort: qsort.cpp
	$(MPI) $(CFLAGS) -o qsort qsort.cpp


clean:
	rm -f ex1 qsort *~ \#* 
