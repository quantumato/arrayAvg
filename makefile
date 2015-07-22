CC = g++
MPICC = mpic++
CFLAGS = -g -Wall -std=c++0x

BIN = arrayAvg
DEP = src/HaloMerged.cpp src/HaloSplit.cpp src/HaloNeighbor.cpp src/HaloNeighbor_nb.cpp

all: $(BIN)

arrayAvg: src/arrayAvg.cpp
	$(MPICC) $(CFLAGS) -o arrayAvg src/arrayAvg.cpp $(DEP)

clean:
	rm -f $(BIN)
