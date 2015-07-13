/*
*
* Program Name: arrayAvg5.cpp
* Author: Kevin Ye
* Date Created: 7/1/17:34
*
*/

#include <iostream>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <stdlib.h>
#include "matrix.h"
#include "Halo6.hpp"

void initialize_matrix(matrix<long int>& A);
void average_local(matrix<long int>& arr1, matrix<long int>& arr2);
//void average_external(long int* arr1, long int* arr2, int n);
void printarray(matrix<long int>& A);
//void setHalo(long int* arr, int rank, int num_p, int length);

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
#ifndef NOTIME
	double start = MPI_Wtime();
#endif
	int lnx, lny; //local number of rows in x and y dimensions
	//NOTE: We are assuming a processor number that is a power of 4
	int gnx = 100, gny = 100; //global number of rows in x and y dimensions

	int rank, size;	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//this is highly dependent on the size of the matrix and the number of processses.
	//I could write this to be more generic but for the purposes of this experiment I'm leaving it like this
	lnx = gnx/sqrt(size) + 2;
	lny = gny/sqrt(size) + 2;

	//matrices to hold values
	matrix<long int> a;
	matrix<long int> b;
	
	a.assign(lnx, lny);
	b.assign(lnx, lny);

	initialize_matrix(a);

	//initialize Halo
	Halo setupHalo(a, rank, size, lnx, lny);
	
	for(int i=0; i<4; i++)
	{
		if(i%2==0)
		{
			setupHalo.Halo_Init(a);
			//average the local-dependent parts of array
			average_local(a, b);	
			/*if(rank == 0)
			{
			std::cout << "A: \n";
			printarray(a);
			std::cout << '\n';
			printarray(b);
			}*/

			//average external-dependent parts of the array
			//average_external(a, b);
			//setupHalo.assignArr(b);

		}
		else
		{
			//if(rank == 0)
				//std::cout << "Starting B\n";
			setupHalo.Halo_Init(b);
			//average the local-depenedent parts of array
			average_local(b, a);
		
			/*if(rank == 0)
			{	
			std::cout << "B: \n";
			printarray(b);
			std::cout << '\n';
			printarray(a);
			}*/
			//setupHalo.Halo_Finalize();
			//average the external-denpendednt parts of the array
			//average_external(a, b, local_len);
			//setupHalo.assignArr(a);
		}
	}

#ifndef NOTIME	
	double finish = MPI_Wtime();
	double time = finish-start;
	double time_sum;
	double time_max;
	double time_min;
	MPI_Allreduce(&time, &time_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&time, &time_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if(rank == 0)
	{
		std::cout << "average= " << time_sum/size << std::endl;
		std::cout << "min= " << time_min << std::endl;
		std::cout << "max= " << time_max << std::endl;
	}
#endif
	
	int the_answer = 42;
	MPI_Status status;
	MPI_Request request;
	if(rank != 0)
		MPI_Irecv(&the_answer, 1, MPI_INT, rank-1, 42, MPI_COMM_WORLD, &request);

	//std::cout << "rank: " << rank << " is printing\n";
	//printarray(a);
	if(rank != 0)
		MPI_Wait(&request, &status);
	if(rank < size-1)
		MPI_Send(&the_answer, 1, MPI_INT, rank+1, 42, MPI_COMM_WORLD);

	MPI_Finalize();

	//std::cout << "time: " << finish-start << std::endl;
	//std::cout <<rank << " The End\n";
	return 0;
}

void average_local(matrix<long int>& arr1, matrix<long int>& arr2)
{
	//assumes n > 1
	for(int i=1; i<arr1.get_nrows()-1; i++)
		for(int j=1; j<arr2.get_ncols()-1; j++)
		//arr2[i][j]=pow(arr1[i-1][j]*arr1[i+1][j]*arr1[i][j+1]*arr1[i][j-1], 0.25);
			arr2[i][j]=(arr1[i-1][j]+arr1[i+1][j]+arr1[i][j+1]+arr1[i][j-1])/4;
}

void initialize_matrix(matrix<long int>& A)
{
	int nx = A.get_nrows();
	int ny = A.get_nrows();
	for(int i=1; i<nx-1; i++)
		for(int j = 1; j<ny-1; j++)
			A[i][j] = rand()%10;
}

/*void average_external(long int* arr1, long int* arr2, int n)
{
	//assumes array minimum size of 3
	arr2[1]=log(exp(arr1[0]+arr1[2]));
	arr2[n-2]=log(exp(arr1[n-3]+arr1[n-1]));
}*/

void printarray(matrix<long int>& A)
{
	int nx = A.get_nrows();
	int ny = A.get_nrows();
	for(int i=0; i<nx; i++)
	{
		for(int j = 0; j<ny; j++)
		std::cout << A[i][j] << " ";
	std::cout << std::endl;
	}
	std::cout << "END" << std::endl;
	return;
}

//******************************************************
//
// Halo Exchange
//
//******************************************************

/*void setupHalo(matrix<long int>& A, int r, int np, int nx, int ny)
{
		int rank; //process ID
		int num_p; //number of processes
		int npx = nx;
		int npy = ny; //number of processors in x and y directions
		int * elementsToRecv; //NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
		int neighbors[4]; //preinitialized to 4 elements TODO: make it dynamically allocated
		int numNeighbors; //number of neighbors
		int totalToSend;
		int totalToReceive;
		MPI_Request* requests;
		MPI_Status *status;

	int local_rows = nx-2;
	int local_cols = ny-2;
	int numToSend = (npx-2)*2 + (npy-2)*2;

	int * elementsToSend = new int[numToSend];
	//receive buffer
	//receive buffer
	requests = new MPI_Request[totalToSend];
													//NOTE: only holds coordinates that change. (i,0), (i,n-1), (0,i), (n-1, i)
	int numNeighbors = 0;
	//calculate the rank of the processor in the processor grid 
	//assume square grid for now
	npx = sqrt(np);
	npy = sqrt(np);

	//calculate the ranks of the neighbors
	//TODO: Rewrite this more efficiently.
	int num_ranks = 0;
	if(rank-1 > 0)
	{
		neighbor[num_ranks] = rank-1;
		num_ranks++;
	}
	if(rank-nx > 0)
	{
		neighbors[num_ranks] = rank-nx;
		num_ranks++;
	}
	if(rank+1<(ny*ny))
	{
		neighbors[num_ranks] = rank+1;
		num_ranks++;
	}
	if(rank+npx<(ny*nx))
	{
		neighbors[num_ranks] = rank+nx;
		num_ranks++;
	}

	//calculate the elements to send	
	//NOTE: Even though this isn't generic this loop will stll hold up if x and y are different dimensions
	int sendIndex = 0;

	for(int i=1; i<nx-1; i++; sendIndex++)
	{
		//send buffer
		elementsToSend[sendIndex] = A[i][1];
		elementsToSend[sendIndex+local_rows] = A[i][local_rows];
	}		
	sendIndex = sendIndex*2 + 2; //NOTE: Possible error source in sendindex being calculated incorrectly

	//TODO: If I make this generic I should include something to keep track of the number of elements		
	for(int i=1; i<ny-1; i++; sendIndex++)
	{
		elementsToSend[sendIndex] = A[1][i];
		elementsToSend[sendIndex+local_rows] = A[local_cols][i];
	}
	//receive buffer
}*/
