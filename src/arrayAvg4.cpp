/*
*
* Program Name: arrayAvg4.cpp
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
#include "Halo4.hpp"
#include <vector>

void initialize_matrix(matrix<long int>& A);
void average_local(matrix<long int>& arr1, matrix<long int>& arr2);
//void average_external(long int* arr1, long int* arr2, int n);
void printarray(matrix<long int>& A);
//void setHalo(long int* arr, int rank, int num_p, int length);

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
#ifndef NOTIME
	std::vector<double> timer;
	double start, finish, diff;
#endif
	int lnx, lny; //local number of rows in x and y dimensions
	//NOTE: We are assuming a processor number that is a power of 4
	//int gnx = 100, gny = 100; //global number of rows in x and y dimensions

	int rank, size;	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//this is highly dependent on the size of the matrix and the number of processses.
	//I could write this to be more generic but for the purposes of this experiment I'm leaving it like this
	if(argc == 3)
	{
	lnx = atoi(argv[1]);
	lny = atoi(argv[2]);
	}
	else
	{
		lnx = 102;
		lny = 102;
	}

	//matrices to hold values
	matrix<long int> a;
	matrix<long int> b;
	
	a.assign(lnx, lny);
	b.assign(lnx, lny);

	initialize_matrix(a);
	initialize_matrix(b);

	//initialize Halo
	Halo setupHalo(a, rank, size, lnx, lny);
	
	for(int i=0; i<2; i++)
	{
		if(i%2==0)
		{
			start = MPI_Wtime();
			setupHalo.Halo_Init(a);
			//average the local-dependent parts of array
			average_local(a, b);	
			finish = MPI_Wtime();
			diff = finish - start;
			//std::cout << diff << std::endl;
			timer.push_back(diff);
		}
		else
		{
			start = MPI_Wtime();
			setupHalo.Halo_Init(b);
			
			//average the local-depenedent parts of array
			average_local(b, a);
			finish = MPI_Wtime();
			diff = finish - start;
			timer.push_back(diff);
		}
	}

#ifndef NOTIME	
	double time;
	double local_min, local_max;
	local_min = timer[0];
	local_max = local_min;
	for(int i=0; i<timer.size(); i++)
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
	double time_min, local_minx, local_maxx;
	MPI_Reduce(&time, &time_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&time, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&time, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	//MPI_Allreduce(&local_min, &local_minx, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	//MPI_Allreduce(&local_max, &local_maxx, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	if(rank == 0)
	{
		std::cout << "average= " << time_sum/size << std::endl;
		std::cout << "min average= " << time_min << std::endl;
		std::cout << "max average= " << time_max << std::endl;
		//std::cout << "min= " << local_minx << std::endl;
		//std::cout << "max= " << local_maxx << std::endl;
	}


#endif
	//MPI_transmission only for the purpose of lining up the output	
	int the_answer = 42;
	MPI_Status status;
	MPI_Request request;
	if(rank != 0)
		MPI_Recv(&the_answer, 1, MPI_INT, rank-1, 42, MPI_COMM_WORLD, &status);

	//std::cout << "rank: " << rank << " is printing\n";
	//printarray(a);
	if(rank < size-1)
		MPI_Send(&the_answer, 1, MPI_INT, rank+1, 42, MPI_COMM_WORLD);

	MPI_Finalize();

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
	for(int i=0; i<nx; i++)
		for(int j = 0; j<ny; j++)
			A[i][j] = rand()%10;
}


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
