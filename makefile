CC = g++
MPICC = mpic++
CFLAGS = -g -Wall -std=c++0x

BIN = arrayAvg
DEP = src/twoDim/HaloMerged.cpp src/twoDim/HaloSplit.cpp src/twoDim/HaloNeighbor.cpp src/twoDim/HaloNeighbor_nb.cpp

all: $(BIN)

arrayAvg: src/twoDim/arrayAvg.cpp
	$(MPICC) $(CFLAGS) -o arrayAvg src/twoDim/arrayAvg.cpp $(DEP)

clean:
	rm -f $(BIN)
