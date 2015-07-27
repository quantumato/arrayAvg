#ifndef HALOMERGED_HPP
#define HALOMERGED_HPP

#include "Halo.hpp"
class HaloMerged : public Halo
{
	private:
		int numNeighbors, recvIndex, sendIndex;
	int neighbors[4]; //max 4 neighbors

		MPI_Request* requests;
		MPI_Status* status;
	public:
		HaloMerged(int r, int np, int nx, int ny);
		~HaloMerged() {};

		bool Halo_Init(matrix<double>& A);
		void Halo_Finalize(matrix<double> &A) {}; //TODO: make inline(?)
};

#endif
