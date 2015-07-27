#include "HaloNeighbor_nb.hpp"
#include<mpi.h>
#include<cmath>
#include<iostream>

HaloNeighbor_nb::HaloNeighbor_nb(int r, int np, int nx, int ny)
{
	//calculate the rank of the processor in the processor grid 
	//assume square grid for now

	//calculate the ranks of the neighbors
	//TODO: Rewrite this mor
	totalSend = local_rows*4;
	totalRecv = 4;

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

bool HaloNeighbor_nb::Halo_Init(matrix<double>& A)
{
	//because local array is a square only one dimension is needed
		//NOTE:getEdgeArray and getExternalArray are not implemented and are purely conceptual
	//for local arrays that are not square use MPI_Neighbor_Ialltoallv
	MPI_Neighbor_Ialltoall(A.getEdgeArray(), local_cols, MPI_DOUBLE, A.getExternalArray(), local_cols, MPI_DOUBLE, cart, &request); 
	return true;
}

void HaloNeighbor_nb::Halo_Finalize()
{
	MPI_Waitall(1, &request, &status);
}

