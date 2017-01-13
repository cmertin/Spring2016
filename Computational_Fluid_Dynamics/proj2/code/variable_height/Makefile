CC = g++
CFLAGS = -O3 -fopenmp --std=c++11 -Wno-unused-result

separation: ground.cpp Dimensions.h
	$(CC) $(CFLAGS) -o ground ground.cpp

clean:
	rm -rf ground *~ \#*
