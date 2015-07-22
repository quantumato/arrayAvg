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
#include <cmath>
#include "matrix.h"

//abstract base class for Halo
class Halo
{
	protected:
		//matrix<long int>& A; //may take up too much memory could just call by reference
		int rank; //process ID
		int num_p; //number of processes
		int npx, npy; //number of processors in x and y directions
		int nix, niy; //number of elements in x and y directions
		int tag;
		//neighbor version needs this
		//long int * elementsToSend; //holds the indices of the elements to be sent
		//NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
		//long int * elementsToRecv; //NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
		int local_rows = npx-2;
		int local_cols = npy-2;
		//merged version requires this
		//int neighbors[4]; //preinitialized to 4 elements TODO: make it dynamically allocated
		int sendLength[4];
		int recvLength[4];
		//merged version requires this
		//int numNeighbors; //number of neighbors
		int totalSend;
		int totalRecv;
		//merged version requires these two
		//MPI_Request* requests;
		//MPI_Status *status;
	public:
		//consolidate initiliazation of common member data
		Halo(int r, int np, int nx, int ny)
		{
			rank = r; //process ID
			num_p = np; //number of processes
			nix = nx;
			niy = ny; 
			npx = sqrt(np);//number of processors in x and y directions. Processors will be a square number
			npy = sqrt(np);
			tag = 42;

			local_rows = niy-2;
			local_cols = nix-2;
			//NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
			//numNeighbors = 0;
			//calculate the rank of the processor in the processor grid 
			//assume square grid for now

			//calculate the ranks of the neighbors
			//TODO: move this into Init
		};

		virtual ~Halo() {}; //virtual destructor required for base class

		virtual bool Halo_Init(matrix<double>& A) =0;
		virtual bool Halo_Finalize() =0;
};
#endif

