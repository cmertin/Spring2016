CC = g++
CFLAGS = -O3 -fopenmp --std=c++11 -Wno-unused-result

train: train.cpp Dimensions.h
	$(CC) $(CFLAGS) -o train train.cpp

clean:
	rm -rf train *~ \#*
