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
		long int* local_array;
		int rank;
		int num_p;
		int length;
		int elementsToSend[2]; //holds the indices of the elements to be sent
		int elementsToRecv[2];
		int neighbors[2];
		int numToSend;
		int total_comm;
		MPI_Request* requests;
		MPI_Status *status;
	public:
		Halo(long int arr[], int r, int np, int len);
		~Halo();

		void Halo_Init();
		void Halo_Finalize();
		void assignArr(long int narray[]);
};
#endif
