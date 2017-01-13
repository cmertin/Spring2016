CC = g++
MPI = mpic++
CFLAGS = -O2 -fopenmp -std=c++11

all: search.cpp bsearch.h ring.cpp
	$(MPI) $(CFLAGS) -o search search.cpp bsearch.h
	$(MPI) $(CFLAGS) -o ring ring.cpp

search: search.cpp bsearch.h
	$(MPI) $(CFLAGS) -o search search.cpp bsearch.h

ring: ring.cpp
	$(MPI) $(CFLAGS) -o ring ring.cpp

clean:
	rm -f ring search *~ \#* 
