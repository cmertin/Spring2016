CC = g++
CFLAGS = -O2 -fopenmp

sort: quicksort.cpp 
	$(CC) $(CFLAGS) -o quicksort quicksort.cpp

scan: genericScan.cpp
	$(CC) $(CFLAGS) -o genericScan genericScan.cpp

all: quicksort.cpp genericScan.cpp
	$(CC) $(CFLAGS) -o quicksort quicksort.cpp
	$(CC) $(CFLAGS) -o genericScan genericScan.cpp

clean:
	rm -f quicksort genericScan *~ \#* 
