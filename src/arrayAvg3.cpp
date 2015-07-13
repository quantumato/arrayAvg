/*
 * Program Name: arrayAvg3.cpp
 *	Author: Kevin Ye
 * Date last modified: 6/25/15
 *
 * This program calculates the averages between two indexes (i-1, i+1) in an integer array and averages
 * them, placing these values in a second array. Then the two arrays are flipped and the proccess is
 * repeated. This version uses MPI and attempts to use Isend to improve performance.
 *
 */

#include <iostream>
#include <algorithm>
#include <cmath>
#include <mpi.h>
#include "Halo3.hpp"
#include <stdlib.h>

void average_local(long int* arr1, long int* arr2, int n);
void average_external(long int* arr1, long int* arr2, int n);
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

	//initialize Halo
	Halo setupHalo(a, rank, size, local_len);
	
	for(int i=0; i<100; i++)
	{
		if(i%2==0)
		{
			setupHalo.Halo_Init();
			//average the local-dependent parts of array
			average_local(a, b, local_len);	

			setupHalo.Halo_Finalize();
			//average external-dependent parts of the array
			average_external(a, b, local_len);
			setupHalo.assignArr(b);

		}
		else
		{
			setupHalo.Halo_Init();
			//average the local-depenedent parts of array
			average_local(b, a, local_len);

			setupHalo.Halo_Finalize();
			//average the external-denpendednt parts of the array
			average_external(a, b, local_len);
			setupHalo.assignArr(a);
		}
	}
	
	double finish = MPI_Wtime();
	double time = finish-start;
	double time_sum;
	MPI_Allreduce(&time, &time_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	if(rank == 0)
		std::cout << "average= " << time_sum/size << std::endl;
	MPI_Finalize();

	std::cout << "time: " << finish-start << std::endl;

	delete[] a;
	delete[] b;
	return 0;
}

void average_local(long int* arr1, long int* arr2, int n)
{
	//assumes n > 1
	for(int i=2; i < n-2; i++)
		arr2[i]=log(exp(arr1[i-1]+arr1[i+1]));
}

void average_external(long int* arr1, long int* arr2, int n)
{
	//assumes array minimum size of 3
	arr2[1]=log(exp(arr1[0]+arr1[2]));
	arr2[n-2]=log(exp(arr1[n-3]+arr1[n-1]));
}

void printarray(long int* arr, int n)
{
	for(int i=0; i<n; i++)
		std::cout << arr[i] << " ";
	std::cout << std::endl;
}

