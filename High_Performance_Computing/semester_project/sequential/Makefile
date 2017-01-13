CC = g++
CFLAGS = -O3 -fopenmp --std=c++11 -Wno-unused-result

sequential: sequential_large.cpp vector.h planet.h hermite.h objects.h
	$(CC) $(CFLAGS) -o sequential sequential_large.cpp

generate: generate.cpp planet.h vector.h
	$(CC) $(CFLAGS) -o generate generate.cpp

clean:
	rm -rf sequential generate *.dat *.pdf *~ \#*
