/*
 *
 * Program Name: arrayAvg.cpp
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
#include "HaloSplit.hpp"
#include "HaloNeighbor.hpp"
#include "HaloMerged.hpp"
#include "HaloNeighbor_nb.hpp"
#include <vector>

void initialize_matrix(matrix<double>& A);
void average_local(matrix<double>& arr1, matrix<double>& arr2);
void printarray(matrix<double>& A);
void average_external(matrix<double>& arr1, matrix<double>& arr2);
bool int_test(char* arg);
//void setHalo(double* arr, int rank, int num_p, int length);

/*
 *
 * params: -s = split
 *					-b = blocking
 *					-n = nonblocking
 *
 */


int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
#ifndef NOTIME
	std::vector<double> timer;
	double start, finish, diff;
#endif
	int lnx, lny; //local number of rows in x and y dimensions
	//NOTE: We are assuming a processor number that is a power of 4

	int rank, size;	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	matrix<double> a;
	matrix<double> b;

	a.assign(lnx, lny);
	b.assign(lnx, lny);

	//this is highly dependent on the size of the matrix and the number of processses.
	//I could write this to be more generic but for the purposes of this experiment I'm leaving it like this
	//TODO: introduce better error checking
	std::string arg1(argv[1]);
	initialize_matrix(a);
	Halo * setupHalo;
	if(argc > 1)
	{
		if(arg1 == "-m")
		{
			setupHalo = new HaloMerged(rank, size, lnx, lny);
		}
		else if(arg1 == "-s")
		{
			setupHalo = new HaloSplit(rank, size, lnx, lny);
		}
		else if(arg1 == "-b") //b for blocking neighbor
		{
			setupHalo = new HaloNeighbor(rank, size, lnx, lny);
		}
		else if(arg1 == "-n") //n for nonblocking neighbor
		{
			setupHalo = new HaloNeighbor_nb(rank, size, lnx, lny);
		}
		else
		{
			std::cout << "Usage: mpirun -np [NUM_PROCESSES] -[HALOTYPE] [x] [y]\n";
			return -1; //no specified reason to return specific number
		}
		for(int i = 2; i<argc; i++)
		{
			if(int_test(argv[i]))
			{
				lnx = atoi(argv[2]);
				lny = atoi(argv[3]);
			}
			else
			{
				std::cout << "Usage: mpirun -np [NUM_PROCESSES] -[HALOTYPE] [x] [y]\n";
				return -2;
			}
		}
	}
	else
	{
		//error message
		std::cerr << "Not enough arguments.\n Usage: mpirun -np [NUM_PROCESSES] -[HALOTYPE] [x] [y]\n";
		return -2;
	}

	//matrices to hold values
	
	initialize_matrix(a);

	matrix<double>* src;
	matrix<double>* dst;

	for(int i=0; i<10; i++)
	{
		//data from src is averaged and loaded into dst then next iteration the arrays are switched.
		if(i%2==0)
		{
			*src = a;
			*dst = b;
		}
		else
		{
			*src = b;
			*dst = a;
		}
		start = MPI_Wtime();
		bool isSplit = setupHalo.Halo_Init(*src);
		if(isSplit)
		{
			//average the local-dependent parts of array
			//TODO: move average functions to matrix?
			average_local(*src, *dst);	
			setupHalo.Halo_Finalize(*src);
		}
		finish = MPI_Wtime();
		diff = finish-start;
		timer.push_back(diff);
		//average external-dependent parts of the array
		average_external(*src, *dst);

	}

#ifndef NOTIME	
	double time;
	double local_min, local_max;
	local_min = timer[0];
	local_max = local_min;

	for(int i = 0; i < timer.size(); i++)
	{
		if(timer[i] < local_min)
			local_min=timer[i];
		if(timer[i] > local_max)
			local_max=timer[i];
		time+=timer[i];	
	}
	time/=timer.size();
	double time_sum;
	double time_max;
	double time_min;
	//double local_minx, local_maxx;
	double time_data[3]; //index 0 = time_sum, index 1 = time_max, index 2 = time_min
	MPI_Op ops[3] = {MPI_SUM, MPI_MIN, MPI_MAX};
	for(int i = 0; i < 3; i++)
	{
		MPI_Reduce(&time, &time_data[i], 1, MPI_DOUBLE, ops[i], 0, MPI_COMM_WORLD);
	}

	//MPI_Reduce(&time, &time_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	//MPI_Reduce(&time, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	//MPI_Reduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	//MPI_Reduce(&local_min, &local_minx, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	//MPI_Reduce(&local_max, &local_maxx, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	if(rank == 0)
	{
		std::cout << "average= " << time_data[0]/size << std::endl;
		std::cout << "max= " << time_max << std::endl;
		std::cout << "min= " << time_min << std::endl;
		//std::cout << "min= " << local_minx << std::endl;
		//std::cout << "max= " << local_maxx << std::endl;
	}

#endif
	//MPI calls to make sure the output lines up to some degree
	//TODO: move to own function
	int the_answer = 42;
	MPI_Status status;
	if(rank != 0)
		MPI_Recv(&the_answer, 1, MPI_INT, rank-1, 42, MPI_COMM_WORLD, &status);

	//std::cout << "rank: " << rank << " is printing\n";
	//printarray(a);
	if(rank < size-1)
		MPI_Send(&the_answer, 1, MPI_INT, rank+1, 42, MPI_COMM_WORLD);

	MPI_Finalize();

	//std::cout << "time: " << finish-start << std::endl;
	//std::cout <<rank << " The End\n";
	return 0;
}

bool int_test(char* arg)
{
	std::string buf = arg;
	for(int i=0; i< buf.size(); i++)
	{
		if(buf[i] < 47 && buf[i] > 58)
			return false;
	}
	return true;
}

void average_local(matrix<double>& arr1, matrix<double>& arr2)
{
	//assumes n > 1
	for(int i=2; i<arr1.get_nrows()-2; i++)
		for(int j=2; j<arr2.get_ncols()-2; j++)
			//arr2[i][j]=pow(arr1[i-1][j]*arr1[i+1][j]*arr1[i][j+1]*arr1[i][j-1], 0.25);
			arr2[i][j]=(arr1[i-1][j]+arr1[i+1][j]+arr1[i][j+1]+arr1[i][j-1])/4;
}

void initialize_matrix(matrix<double>& A)
{
	int nx = A.get_nrows();
	int ny = A.get_ncols();
	for(int i=0; i<nx; i++)
		for(int j = 0; j<ny; j++)
			A[i][j] = rand()%10;
}

void average_external(matrix<double>& arr1, matrix<double>& arr2)
{
	//assumes array minimum size of 3
	int nx = arr1.get_nrows();
	int ny = arr1.get_ncols();
	for(int i=1; i<nx-1; i++)
	{
		arr2[1][i]=(arr1[0][i]+arr1[2][i]+arr1[1][i-1]+arr1[1][i+1])/4;
		arr2[nx-2][i]=(arr1[nx-3][i]+arr1[nx-1][i]+arr1[nx-2][i-1]+arr1[nx-2][i+1])/4;
	}
	for(int i = 1; i<ny-1; i++)
	{
		arr2[i][1]=(arr1[i][0]+arr1[i][2]+arr1[i-1][1]+arr1[i+1][1])/4;
		arr2[i][nx-2]=(arr1[i][nx-3]+arr1[i][nx-1]+arr1[i+1][nx-2]+arr1[i-1][nx-2])/4;
	}

}

void printarray(matrix<double>& A)
{
	int nx = A.get_nrows();
	int ny = A.get_nrows();
	for(int i=0; i<nx; i++)
	{
		for(int j = 0; j<ny; j++)
			std::cout << A[i][j] << " ";
		std::cout << std::endl;
	}
	std::cout << '\n';
	return;
}

//******************************************************
//
// Halo Exchange
//
//******************************************************

/*void setupHalo(matrix<double>& A, int r, int np, int nx, int ny)
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
