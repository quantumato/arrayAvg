#ifndef HALOMERGED_HPP
#define HALOMERGED_HPP

#include "Halo.hpp"
class HaloMerged : public Halo
{
	private:
		long int * elementsToSend; //holds the indices of the elements to be sent
													//NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
		long int * elementsToRecv; //NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
		int numNeighbors, recvIndex, sendIndex;
	int neighbors[4]; //max 4 neighbors

		MPI_Request* requests;
	public:
		HaloMerged(int r, int np, int nx, int ny);
		~HaloMerged() {};

		bool Halo_Init(matrix<double>& A);
		void Halo_Finalize(matrix<double> &A) {}; //TODO: make inline(?)
};

#endif
