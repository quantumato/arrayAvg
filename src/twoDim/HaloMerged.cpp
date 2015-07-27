#include "HaloMerged.hpp"
#include "Halo.hpp"
#include<mpi.h>
#include<cmath>
#include<iostream>

HaloMerged::HaloMerged(int r, int np, int nx, int ny)
	: Halo{r, np, nx, ny}
{
	rank = r; //process ID
	num_p = np; //number of processes
	nix = nx;
	niy = ny; 
	npx = sqrt(np);//number of processors in x and y directions. Processors will be a square number
	npy = sqrt(np);
	tag = 42;

	local_rows = ny-2;
	local_cols = nx-2;
	numNeighbors = 0; //used as 

	//calculate the rank of the processor in the processor grid 
	//assume square grid for now

	//calculate the ranks of the neighbors
	//TODO: move this into Init

	totalSend=0;
	recvIndex = 0;
	if(rank >= npx) //check if not on top
	{
		neighbors[numNeighbors] = rank-npx;
		numNeighbors++;
		totalSend+=local_cols;
		//if local array is not a square send/recvLength would vary between these ifs
		sendLength[recvIndex]=local_cols;
		recvLength[recvIndex]=local_cols;
		recvIndex++;

	}
	if(rank%npx > 0) //check if not on left side
	{
		neighbors[numNeighbors] = rank-1;
		numNeighbors++;
		totalSend+=local_rows;
		sendLength[recvIndex]=local_cols;
		recvLength[recvIndex]=local_cols;
		recvIndex++;
	}
	if(rank+npx<(npx*npy)) //check not on bottom side
	{
		neighbors[numNeighbors] = rank+npx;
		numNeighbors++;
		totalSend+=local_cols;
		sendLength[recvIndex]=local_cols;
		recvLength[recvIndex]=local_cols;
		recvIndex++;

	}
	if(rank%npx!=npx-1) //check not on right side
	{
		neighbors[numNeighbors] = rank+1;
		numNeighbors++;
		totalSend+=local_rows;
		sendLength[recvIndex]=local_cols;
		recvLength[recvIndex]=local_cols;
		recvIndex++;

	}
}

bool HaloMerged::Halo_Init(matrix<double>& A)
{
		requests = new MPI_Request[numNeighbors];
	status = new MPI_Status[numNeighbors];

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
		MPI_Send(A.getExternal(i), n_send, MPI_DOUBLE, neighbors[i], tag, MPI_COMM_WORLD);
	}
	MPI_Waitall(numNeighbors, requests, status);

	delete [] status;
	delete [] requests;
	return false;
}
