/*
* Program Name: arrayAvg1.cpp
*	Author: Kevin Ye
* Date last modified: 6/23/15
*
* this program calculates the averages between two indexes (i-1, i+1) in an integer array and averages
* them, placing these values in a second array. Then the two arrays are flipped and the proccess is
* repeated.
*
*/


#include <iostream>
#include <algorithm>
#include <cmath>

void average(int* arr1, int* arr2, int n);
void printarray(int* arr, int n);
int main(int argc, char* argv[])
{
	int size = 1000000;
	int * a = new int[size];
	int * b = new int[size];
	
	for(int i=0; i<size; i++)
		a[i]=2*(i+1);
	for(int i=0; i < 100; i++)
	{
		if(i%2==0)
			average(a, b, size);	
		else
			average(b, a, size);
	}
	//printarray(b, size);
	return 0;
}

void average(int* arr1, int* arr2, int n)
{
	//assumes n > 1
	arr2[0] = (arr1[1] + 1)/2;
	arr2[n-1]= (arr1[n-2] + 42)/2;
	for(int i=1; i < n-1; i++)
		arr2[i]=sqrt((arr1[i-1]*arr1[i+1])/2);
}

void printarray(int* arr, int n)
{
	for(int i=0; i<n; i++)
		std::cout << arr[i] << " ";
	std::cout << std::endl;
}

