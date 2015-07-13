CC = g++
MPICC = mpic++
CFLAGS = -g -Wall -std=c++0x

BIN = arrayAvg1 arrayAvg2 arrayAvg3 arrayAvg4 arrayAvg5 arrayAvg6

all: $(BIN)

arrayAvg1: src/arrayAvg1.cpp
	$(CC) $(CFLAGS) -o arrayAvg1 src/arrayAvg1.cpp

arrayAvg2: src/arrayAvg2.cpp
	$(MPICC) $(CFLAGS) -o arrayAvg2 src/arrayAvg2.cpp

arrayAvg3: src/arrayAvg3.cpp
	$(MPICC) $(CFLAGS) -o arrayAvg3 src/arrayAvg3.cpp src/Halo3.cpp

arrayAvg4: src/arrayAvg4.cpp src/Halo4.cpp
	$(MPICC) $(CFLAGS) -o arrayAvg4 src/arrayAvg4.cpp src/Halo4.cpp

arrayAvg5: src/arrayAvg5.cpp src/Halo5.cpp
	$(MPICC) $(CFLAGS) -o arrayAvg5 src/arrayAvg5.cpp src/Halo5.cpp

arrayAvg6: src/arrayAvg6.cpp src/Halo6.cpp
	$(MPICC) $(CFLAGS) -o arrayAvg6 src/arrayAvg6.cpp src/Halo6.cpp

clean:
	rm -f $(BIN)
