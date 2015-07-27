/*
*
* Program Name: setupHalo.hpp
* Author: Kevin Ye
* Date Last Modifed: 6/25/15
*
*/

#ifndef HALO3_HPP
#define HALO3_HPP
#include <mpi.h>
#include "matrix.h"

class Halo
{
	private:
		//matrix<long int>& A; //may take up too much memory could just call by reference
		int rank; //process ID
		int num_p; //number of processes
		int npx, npy; //number of processors in x and y directions
		int nix, niy; //number of elements in x and y directions
		int tag;
		long int * elementsToSend; //holds the indices of the elements to be sent
													//NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
		long int * elementsToRecv; //NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
		int local_rows = npx-2;
		int local_cols = npy-2;
		int neighbors[4]; //preinitialized to 4 elements TODO: make it dynamically allocated
		int sendLength[4];
		int recvLength[4];
		int numNeighbors; //number of neighbors
		int totalSend;
		int totalRecv;
		MPI_Request* requests;
		MPI_Status *status;
	public:
		Halo(matrix<long int>& nA, int r, int np, int nx, int ny);
		~Halo();

		void Halo_Init(matrix<long int>& A);
		void Halo_Finalize(matrix<long int>& A);
		//void assignArr(matrix<long int>& nA);
};
#endif

