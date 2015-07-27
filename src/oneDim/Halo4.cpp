#include "Halo4.hpp"
#include<mpi.h>
#include<cmath>
#include<iostream>

Halo::Halo(matrix<long int>& A, int r, int np, int nx, int ny)
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
	//NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
	numNeighbors = 0;
	//calculate the rank of the processor in the processor grid 
	//assume square grid for now

	//calculate the ranks of the neighbors
	//TODO: Rewrite this more efficiently.

	totalSend=0;
	if(rank >= npx) //check if not on top
	{
		neighbors[numNeighbors] = rank-npx;
		numNeighbors++;
		totalSend+=local_cols;
	}
	if(rank%npx > 0) //check if not on left side
	{
		neighbors[numNeighbors] = rank-1;
		numNeighbors++;
		totalSend+=local_rows;
	}
	if(rank+npx<(npx*npy)) //check not on bottom side
	{
		neighbors[numNeighbors] = rank+npx;
		numNeighbors++;
		totalSend+=local_cols;
	}
	if(rank%npx!=npx-1) //check not on right side
	{
		neighbors[numNeighbors] = rank+1;
		numNeighbors++;
		totalSend+=local_rows;
	}
	totalRecv = 0; //REDUNDANT see numNeighbors for function

	elementsToSend = new long int[totalSend];
	elementsToRecv = new long int[totalSend];
	//calculate the elements to send	
	//NOTE: Even though this isn't generic this loop will stll hold up if x and y are different dimensions

	//If statements to check on 
	//TODO: improve this if structure if possible
	//checking if rank is on edge of the square of ranks
	int sendIndex = 0;
	if(rank >= npx) //check if not on top edge
	{
		for(int i=1; i<nix-1; i++, sendIndex++)
		{
			//send buffer
			elementsToSend[sendIndex] = A[1][i];
		}		
		//totalSend scales with number of processes
		sendLength[totalRecv]=local_cols;
		recvLength[totalRecv]=local_cols;
		totalRecv++;
	}

	if(rank%npx!=0) //check if not on left edge
	{
		for(int i=1; i<nix-1; i++, sendIndex++)
		{
			elementsToSend[sendIndex] = A[i][1];
		}			
		sendLength[totalRecv]=local_rows;
		recvLength[totalRecv]=local_rows;
		totalRecv++;
	}

	if(rank< num_p-npx) //check if not on bottom edge
	{
		for(int i=1; i<niy-1; i++, sendIndex++)
		{
			elementsToSend[sendIndex] = A[local_rows][i];
		}
		sendLength[totalRecv]=local_cols;
		recvLength[totalRecv]=local_cols;
		totalRecv++;
	}

	if(rank%npx!=npx-1) //check if not on right edge
	{
		for(int i=1; i<niy-1; i++, sendIndex++)
		{
			elementsToSend[sendIndex] = A[i][local_rows];
		}
		sendLength[totalRecv]=local_rows;
		recvLength[totalRecv]=local_rows;
		totalRecv++;
	}
}

//
//
// Halo Exchange
//
//

void Halo::Halo_Init(matrix<long int>& A)
{
	//post Irecv	

	//TODO: replace totalRecv with numNeighbors
	requests = new MPI_Request[totalRecv];
	status = new MPI_Status[totalRecv];

	long int * n = (long int *) elementsToRecv;

	for(int i=0; i<totalRecv; i++)
	{
		int n_recv = recvLength[i];
		MPI_Irecv(n, n_recv, MPI_LONG, neighbors[i], tag, MPI_COMM_WORLD, requests+i); 
		
		n+=n_recv;
	}
	//Sends
	long int * m = (long int *) elementsToSend;
	for(int i=0; i<totalRecv; i++)
	{
		int n_send = sendLength[i];
		MPI_Send(m, n_send, MPI_LONG, neighbors[i], tag, MPI_COMM_WORLD);
		m+=n_send;
	}
	MPI_Waitall(totalRecv, requests, status);
	
	//move received elements to array
	int recvIndex = 0;
	if(rank >= npx) //check if on top edge
	{
		for(int i=1; i<nix-1; i++, recvIndex++)
		{
			//send buffer
			A[0][i] = elementsToRecv[recvIndex];
		}		
		//totalSend scales with number of processes
	}

	if(rank%npx!=0) //check if on left edge
	{
		for(int i=1; i<nix-1; i++, recvIndex++)
		{
			A[i][0] = elementsToRecv[recvIndex];
		}			
	}

	if(rank< num_p-npx) //check if on bottom edge
	{
		for(int i=1; i<niy-1; i++, recvIndex++)
		{
			A[local_rows+1][i] = elementsToRecv[recvIndex];
		}
	}

	if(rank%npx!=npx-1) //check if on right edge
	{
		for(int i=1; i<niy-1; i++, recvIndex++)
		{
			A[i][local_cols+1] = elementsToRecv[recvIndex];
		}
	}
	delete [] status;
	delete [] requests;
}

Halo::~Halo()
{
	delete [] elementsToSend;
	delete [] elementsToRecv;
}

/*void Halo::assignArr(matrix<long int>& nA)
{
	A = na;
}*/
