/*
*
* Program Name: Halo3.cpp
* Author: Kevin Ye
* Date Last Modified: 6/26/15
* 
* Implements Halo Exchange for arrayAvg3.cpp
*
*/
#include "Halo3.hpp"

Halo::Halo(long int arr[], int r, int np, int len)
{
	local_array = arr;	
	rank = r;
	num_p = np;
	length = len;
	
	numToSend = 0;

	if(rank > 0)
	{
		elementsToSend[numToSend] = 1;
		elementsToRecv[numToSend] = 0;
		neighbors[numToSend] = rank-1; //keeps track of what process this data is sent to
		numToSend++;
	}
	else
		arr[0]=arr[1];
	if(rank < num_p-1)
	{
		elementsToSend[numToSend] = length-2;
		elementsToRecv[numToSend] = length-1;
		neighbors[numToSend] = rank+1;
		numToSend++;
	}
	else
		arr[length-1] = arr[length-2];
	//total number of communications
	total_comm = numToSend*2;
	requests = new MPI_Request[total_comm];
}

Halo::~Halo()
{
	delete [] requests;
	delete [] status;
}

void Halo::Halo_Init()
{
	int tag = 42;
	for(int i=0; i < numToSend; i++)
		MPI_Irecv(local_array + elementsToRecv[i], 1, MPI_LONG, neighbors[i], tag, MPI_COMM_WORLD, requests+i);
	for(int i=0; i < numToSend; i++)
		MPI_Isend(local_array + elementsToSend[i], 1, MPI_LONG, neighbors[i], tag, MPI_COMM_WORLD, requests+i+numToSend);
}

void Halo::Halo_Finalize()
{
		status = new MPI_Status[total_comm];
		MPI_Waitall(total_comm, requests, status);
}

void Halo::assignArr(long int narray[])
{
	local_array = narray;
	if(rank == 0)
		local_array[0]=local_array[1];
	if(rank == num_p-1)
		local_array[length-1]=local_array[length-2];
}
