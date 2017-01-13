CC = g++
CFLAGS = -O3 --std=c++11 -Wno-unused-result

3heat: 3heat.cpp
	$(CC) $(CFLAGS) -o 3heat 3heat.cpp

2c: 2c.cpp
	$(CC) $(CFLAGS) -o 2c 2c.cpp

2b: 2b.cpp
	$(CC) $(CFLAGS) -o 2b 2b.cpp

all: 3heat.cpp 2b.cpp 2c.cpp
	$(CC) $(CFLAGS) -o 3heat 3heat.cpp
	$(CC) $(CFLAGS) -o 2b 2b.cpp
	$(CC) $(CFLAGS) -o 2c 2c.cpp

clean:
	rm -rf 2b 2c 3heat *.dat *.pdf *~ \#*
