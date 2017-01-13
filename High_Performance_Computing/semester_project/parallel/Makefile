CC = g++
MPICC = mpic++
CFLAGS = -O3 -fopenmp --std=c++11 -Wno-unused-result

parallel: src/barneshut.cpp src/vector.h src/planet.h src/hermite.h 
	$(MPICC) $(CFLAGS) -o barneshut src/barneshut.cpp

generate: src/generate.cpp src/planet.h src/vector.h
	$(CC) $(CFLAGS) -o generate src/generate.cpp

clean:
	rm -rf barneshut generate *.pdf *~ \#* src/*~ src/\#* #data/*.dat data/*.pdf
