#include "HaloSplit.hpp"
#include<mpi.h>
#include<cmath>
#include<iostream>

HaloSplit::HaloSplit(int r, int np, int nx, int ny)
{
	//calculate the rank of the processor in the processor grid 
	//assume square grid for now

	//calculate the ranks of the neighbors
	//TODO: Rewrite this more efficiently.

	totalSend=0;
	numNeighbors = 0;
	if(rank >= npx) //check if not on top
	{
		neighbors[numNeighbors] = rank-npx;
		numNeighbors++;
		totalSend+=local_cols;
		sendLength[recvIndex]=local_cols;
		recvLength[recvIndex]=local_cols;

	}
	if(rank%npx > 0) //check if not on left side
	{
		neighbors[numNeighbors] = rank-1;
		numNeighbors++;
		totalSend+=local_rows;
		sendLength[recvIndex]=local_rows;
		recvLength[recvIndex]=local_rows;

	}
	if(rank+npx<(npx*npy)) //check not on bottom side
	{
		neighbors[numNeighbors] = rank+npx;
		numNeighbors++;
		totalSend+=local_cols;
		sendLength[recvIndex]=local_rows;
		recvLength[recvIndex]=local_rows;

	}
	if(rank%npx<npx-1) //check not on right side
	{
		neighbors[numNeighbors] = rank+1;
		numNeighbors++;
		totalSend+=local_rows;
		sendLength[recvIndex]=local_rows;
		recvLength[recvIndex]=local_rows;

	}
	numNeighbors = 0;
}

bool HaloSplit::Halo_Init(matrix<double>& A)
{
	//post Irecv	

	//TODO: replace numNeighbors with totalRecv
	requests = new MPI_Request[2*numNeighbors];
	status = new MPI_Status[2*numNeighbors];

	for(int i=0; i<numNeighbors; i++)
	{
		int n_recv = recvLength[i];
		//made up some syntax to deal with edges. getEdge is implemented for conceptual purposes
		//and will not actually work
		MPI_Irecv(A.getEdge(i), n_recv, MPI_DOUBLE, neighbors[i], tag, MPI_COMM_WORLD, requests+i); 
	}
	//Sends
	for(int i=0; i<numNeighbors; i++)
	{
		int n_send = sendLength[i];
		MPI_Isend(A.getExternals(i), n_send, MPI_DOUBLE, neighbors[i], tag, MPI_COMM_WORLD, requests+i+numNeighbors);
	}
	return true;
}

//separate wait function
void HaloSplit::Halo_Finalize(matrix<double int>& A)
{
	MPI_Waitall(numNeighbors, requests, status);
	delete [] status;
	delete [] requests;
}

/*void Halo::assignArr(matrix<long int>& nA)
{
	A = na;
}*/
