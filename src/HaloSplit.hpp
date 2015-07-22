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
		//matrix<long int>& A; //may take up too much memory could just call by reference
		//long int * elementsToSend; //holds the indices of the elements to be sent
													//NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
		//long int * elementsToRecv; //NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
		int neighbors[4];
		int sendLength[4];
		int recvLength[4];
		int totalSend;
		int totalRecv;
		MPI_Request* requests;
		MPI_Status *status;
	public:
		HaloSplit(matrix<double>& nA, int r, int np, int nx, int ny);
		~HaloSplit() {};

		bool Halo_Init(matrix<double>& A);
		void Halo_Finalize(matrix<double>& A);
		//void assignArr(matrix<long int>& nA);
};
#endif

