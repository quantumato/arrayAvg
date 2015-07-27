#include "HaloNeighbor.hpp"
#include<mpi.h>
#include<cmath>
#include<iostream>

HaloNeighbor::HaloNeighbor(int r, int np, int nx, int ny)
	: Halo(r, np, nx, ny)
{
	//assume square local matrix
	//calculate the ranks of the neighbors
	//TODO: Rewrite this mor
	totalSend = local_rows*4;
	numNeighbors = 0; //REDUNDANT see numNeighbors for function


	//calculate the elements to send	
	//NOTE: Even though this isn't generic this loop will stll hold up if x and y are different dimensions

	//If statements to check on 
	//pull all the elements that are to be sent into one array

	//std::cout << "calculated elements to Send\n";
	//TODO: move these somewhere more relevant

	//TODO:move these to the header file
	ndim[0] = npy;
	ndim[1] = npx;
	int period[2];
	period[0]=false;
	period[1]=false;
	int reorder = false;

	//
	// Create cartesian coordinate grid
	//
	MPI_Cart_create(MPI_COMM_WORLD, 2, ndim, period, reorder, &cart);
	int coord[2];
	MPI_Cart_coords(cart, rank, 2, coord);
}

bool HaloNeighbor::Halo_Init(matrix<double>& A)
{
	//NOTE: getEdgeArray and getExternalArray are not implemented and are purely conceptual
	//for a local array that is not square use MPI_Neighbor_alltoallv
	MPI_Neighbor_alltoall(A.getEdgeArray(), local_cols, MPI_DOUBLE, A.getExternalArray(), local_cols, MPI_LONG, cart); 
	return false;
}

