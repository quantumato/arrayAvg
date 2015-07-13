/*
* Program Name: arrayAvg2.cpp
*	Author: Kevin Ye
* Date last modified: 6/25/15
*
* This program calculates the averages between two indexes (i-1, i+1) in an integer array and averages
* them, placing these values in a second array. Then the two arrays are flipped and the proccess is
* repeated. This version uses MPI.
*
*/

#include <iostream>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include <stdlib.h>

void average(long int* arr1, long int* arr2, int n);
void printarray(long int* arr, int n);
void setHalo(long int* arr, int rank, int num_p, int length);

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	double start = MPI_Wtime();
	int local_len;
	int global_len = 100000000;

	int rank, size;	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == size-1)
		local_len = global_len/size + global_len%size + 2;
	else
		local_len = global_len/size+2;
	//external values are stored in index 0 and index n-1. The local data is 1 -> n-2

	long int * a = new long int[local_len];
	long int * b = new long int[local_len];
	//initialize array
	for(int i=1; i<local_len-1; i++)
		a[i]=rand()%10;

	for(int i=0; i<100; i++)
	{
		if(i%2==0)
		{
			//set up the halo of the array
			setHalo(a, rank, size, local_len);

			//average the array
			average(a, b, local_len);	
			//std::cout << "process: " << rank << " b: ";
			//printarray(b, local_len);
		}
		else
		{
			setHalo(b, rank, size, local_len);
			average(b, a, local_len);
			//std::cout << "process: " << rank << " a: ";
			//printarray(a, local_len);
		}
	}

	double finish = MPI_Wtime();
	double time = finish-start;
	double time_sum;
	MPI_Allreduce(&time, &time_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if(rank == 0)
		std::cout << "average: " << time_sum/size << std::endl;
	std::cout << "time: " << finish-start << std::endl;
	MPI_Finalize();
	delete[] a;
	delete[] b;
	return 0;
}

void average(long int* arr1, long int* arr2, int n)
{
	//assumes n > 1
	for(int i=1; i < n-1; i++)
		arr2[i]=log(exp(arr1[i-1]+arr1[i+1]));
}

void printarray(long int* arr, int n)
{
	for(int i=0; i<n; i++)
		std::cout << arr[i] << " ";
	std::cout << std::endl;
}

void setHalo(long int* arr, int rank, int num_p, int length)
{
	int tag = 42;
	MPI_Status status;
	//check for not being the first array
	if(rank > 0)
	{
		int nrank = rank-1;
		MPI_Request request;

		MPI_Irecv(arr, 1, MPI_LONG, nrank, tag, MPI_COMM_WORLD, &request);
		MPI_Send(arr+1, 1, MPI_LONG, nrank, tag, MPI_COMM_WORLD);
		MPI_Wait(&request, &status);
	}
	else
		arr[0]=arr[1];

	if(rank < (num_p-1))
	{
		int nrank = rank + 1;
		MPI_Request request;

		MPI_Irecv(arr+length-1, 1, MPI_LONG, nrank, tag, MPI_COMM_WORLD, &request);
		MPI_Send(arr+length-2, 1, MPI_LONG, nrank, tag, MPI_COMM_WORLD);
		MPI_Wait(&request, &status);
	}
	else
		arr[length-1]=arr[length-2];
}
