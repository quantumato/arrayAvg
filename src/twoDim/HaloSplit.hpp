/*
*
* Program Name: setupHalo.hpp
* Author: Kevin Ye
* Date Last Modifed: 6/25/15
*
*/

#ifndef HALOSPLIT_HPP
#define HALOSPLIT_HPP
#include <mpi.h>
#include "matrix.h"
#include "Halo.hpp"

class HaloSplit : public Halo
{
	private:
		int numNeighbors;
		int neighbors[4];
		int sendLength[4];
		int recvLength[4];
		int totalSent, totalRecv, recvIndex;
		MPI_Request* requests;
		MPI_Status *status;
	public:
		HaloSplit(int r, int np, int nx, int ny);
		~HaloSplit() {};

		bool Halo_Init(matrix<double>& A);
		void Halo_Finalize(matrix<double>& A);
};
#endif

