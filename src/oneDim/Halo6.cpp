#include "Halo6.hpp"
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
	//numNeighbors; //number of neighbors
	tag = 42;

	local_rows = ny-2;
	local_cols = nx-2;
	numNeighbors = 0;
	//calculate the rank of the processor in the processor grid 
	//assume square grid for now

	//calculate the ranks of the neighbors
	//TODO: Rewrite this mor
	totalSend = local_rows*4;
	totalRecv = 0; //REDUNDANT see numNeighbors for function

	elementsToSend = new long int[totalSend];
	elementsToRecv = new long int[totalSend];

	//calculate the elements to send	
	//NOTE: Even though this isn't generic this loop will stll hold up if x and y are different dimensions

	//If statements to check on 
	//pull all the elements that are to be sent into one array
	int sendIndex = 0;

		for(int i=1; i<nix-1; i++, sendIndex++)
		{
			elementsToSend[sendIndex] = A[i][1];
		}			
		sendLength[totalRecv]=local_rows;
		recvLength[totalRecv]=local_rows;
		totalRecv++;

		for(int i=1; i<niy-1; i++, sendIndex++)
		{
			elementsToSend[sendIndex] = A[i][local_rows];
		}
		sendLength[totalRecv]=local_rows;
		recvLength[totalRecv]=local_rows;
		totalRecv++;

		for(int i=1; i<nix-1; i++, sendIndex++)
		{
			//send buffer
			elementsToSend[sendIndex] = A[1][i];
		}		
		//totalSend scales with number of processes
		sendLength[totalRecv]=local_cols;
		recvLength[totalRecv]=local_cols;
		totalRecv++;

		for(int i=1; i<niy-1; i++, sendIndex++)
		{
			elementsToSend[sendIndex] = A[local_rows][i];
		}
		sendLength[totalRecv]=local_cols;
		recvLength[totalRecv]=local_cols;
		totalRecv++;
	sdisp[0] = sendLength[0];
	rdisp[0] = recvLength[0];
	//set up arrays that hold displacement info
	for(int i=1; i < totalRecv; i++)
	{
		sdisp[i]=sendLength[i]+sdisp[i-1];
		rdisp[i]=sendLength[i]+rdisp[i-1];
	}

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
	MPI_Neighbor_alltoall(elementsToSend, sendLength[0], MPI_LONG, elementsToRecv, recvLength[0], MPI_LONG, cart); 
	int recvIndex = 0;
		for(int i=1; i<nix-1; i++, recvIndex++)
		{
			//send buffer
			A[0][i] = elementsToRecv[recvIndex];
			//A[0][i] = 96;
		}		
		//totalSend scales with number of processes

		for(int i=1; i<nix-1; i++, recvIndex++)
		{
			A[i][0] = elementsToRecv[recvIndex];
			//A[i][0] = 69;
		}			

		for(int i=1; i<niy-1; i++, recvIndex++)
		{
			A[local_rows+1][i] = elementsToRecv[recvIndex];
			//A[local_rows+1][i] = 42;
		}

		for(int i=1; i<niy-1; i++, recvIndex++)
		{
			A[i][local_cols+1] = elementsToRecv[recvIndex];
			//A[i][local_rows+1] = 21;
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
