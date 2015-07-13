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
	//neighbors[4]; //preinitialized to 4 elements TODO: make it dynamically allocated
	//numNeighbors; //number of neighbors
	tag = 42;

	//std::cout << "npx=" << npx << std::endl;
	//std::cout << "npy=" << npy << std::endl;

	local_rows = ny-2;
	local_cols = nx-2;
	//NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
	numNeighbors = 0;
	//calculate the rank of the processor in the processor grid 
	//assume square grid for now

	//calculate the ranks of the neighbors
	//TODO: Rewrite this more efficiently.

	totalSend=0;
	if(rank%npx > 0) //check if not on left side
	{
		neighbors[numNeighbors] = rank-1;
		numNeighbors++;
		totalSend+=local_rows;
	}
	if(rank%npx<npx-1) //check not on right side
	{
		neighbors[numNeighbors] = rank+1;
		numNeighbors++;
		totalSend+=local_rows;
	}
	if(rank >= npx) //check if not on top
	{
		neighbors[numNeighbors] = rank-npx;
		numNeighbors++;
		totalSend+=local_cols;
	}
	if(rank+npx<(npx*npy)) //check not on bottom side
	{
		neighbors[numNeighbors] = rank+npx;
		numNeighbors++;
		totalSend+=local_cols;
	}

	totalRecv = 0; //REDUNDANT see numNeighbors for function

	elementsToSend = new long int[100];
	elementsToRecv = new long int[100];

//*****************************************************
//TEMPORARY DEBUG CODE
//Initialize elementsToSend/elementsToRecv
	for(int i=0; i<100; i++)
	{
		elementsToSend[i] = 0;
		elementsToRecv[i] = 0;
	}
//END DEBUG CODE
//*****************************************************

	//calculate the elements to send	
	//NOTE: Even though this isn't generic this loop will stll hold up if x and y are different dimensions

	//If statements to check on 
	//TODO: improve this if structure if possible
	//checking if rank is on edge of the square of ranks
	int sendIndex = 0;
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
	sdisp[0] = sendLength[0];
	rdisp[0] = recvLength[0];
	//set up arrays that hold displacement info
	for(int i=1; i < totalRecv; i++)
	{
		sdisp[i]=sendLength[i]+sdisp[i-1];
		rdisp[i]=sendLength[i]+rdisp[i-1];
	}
	//std::cout << "totalRecv: " << totalRecv << "\n";
	//for(int i = 0; i < totalRecv; i++)
	//std::cout << recvLength[i] << " ";
	//std::cout << '\n';

	/*for(int i=0; i<totalSend; i++)
		std::cout << elementsToSend[i] << " ";
		std::cout << std::endl;*/

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
	std::cout << "rank: " << rank << " coords: " << coord[0] << " " << coord[1] << '\n';
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



	//std::cout << "rank: " << rank << " ";
	//for(int i=0; i < totalRecv; i++)
	//std::cout << neighbors[i] << " ";
	//std::cout << '\n';

	/*	std::cout << "totalRecv: " << totalRecv << '\n';
			std::cout << "RecvList: ";
			for(int i=0; i<totalSend; i++)
			{
			std::cout << elementsToRecv[i] << " ";
			}
			std::cout << std::endl;*/

	/*for(int i=0; i<totalRecv; i++)
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
	std::cout << "SENDING: ";
	for(int i=0; i<totalSend; i++)
	std::cout << *(m+i) << " ";
	std::cout << std::endl;
	}
	for(int i=0; i<totalRecv; i++)
	{
	MPI_Wait(requests+i, status);
	std::cout << "RECEIVING: ";
	for(int i=0; i<totalSend; i++)
	std::cout << elementsToRecv[i] << " ";
	std::cout << " neighbor: " << neighbors[i];
	std::cout << '\n';*/
	/*if(rank == 0)
		{
		std::cout << "wait from " << neighbors[i] << " complete.\n";
		}
		}	*/
	/*std::cout << "rank: " << rank << " Elements sent: ";
		for(int i=0; i<totalSend; i++)
		std::cout << elementsToSend[i] << " ";
		std::cout << std::endl;*/
	//MPI_Waitall(totalRecv, requests, status);
	//move received elements to array
	//TODO: displacement array
	/*std::cout << "rank: " << rank << " ";
		for(int i=0; i < totalSend; i++)
		std::cout << elementsToRecv[i] << " ";
		std::cout << '\n';*/
	MPI_Neighbor_alltoall(elementsToSend, sendLength[0], MPI_LONG, elementsToRecv, recvLength[0], MPI_LONG, cart); 
	std::cout << "rank: " << rank << " RECEIVING: ";
	for(int i=0; i<100; i++)
	std::cout << elementsToRecv[i] << " ";
	std::cout << '\n';
	int recvIndex = 0;
	if(rank >= npx) //check if on top edge
	{
		for(int i=1; i<nix-1; i++, recvIndex++)
		{
			//send buffer
			A[0][i] = elementsToRecv[recvIndex];
			//A[0][i] = 96;
		}		
		//totalSend scales with number of processes
	}

	if(rank%npx!=0) //check if on left edge
	{
		for(int i=1; i<nix-1; i++, recvIndex++)
		{
			A[i][0] = elementsToRecv[recvIndex];
			//A[i][0] = 69;
		}			
	}

	if(rank< num_p-npx) //check if on bottom edge
	{
		for(int i=1; i<niy-1; i++, recvIndex++)
		{
			A[local_rows+1][i] = elementsToRecv[recvIndex];
			//A[local_rows+1][i] = 42;
		}
	}

	if(rank%npx!=npx-1) //check if on right edge
	{
		for(int i=1; i<niy-1; i++, recvIndex++)
		{
			A[i][local_cols+1] = elementsToRecv[recvIndex];
			//A[i][local_rows+1] = 21;
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
