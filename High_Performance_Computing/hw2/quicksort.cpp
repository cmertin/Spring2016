/*
  Author: Christopher Mertin
  Date: January 16, 2016
  File: quicksort.cpp
  This program first builds an array of random points, and then sorts them using
  quicksort. Each point has an x and y value, and this program first sorts on
  the y values, and if they're equal, then it sorts on the x values. It is also
  generalized so that it works with integers, longs, floats, and doubles as well.
*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cmath>
#include <omp.h>
#include <algorithm>
#include "points.h"

using namespace std;


template <typename T>
void QuickSort(T *A, int num, size_t size, int (*Comparison)(T,T), bool parallel=true);
template <typename T>
void QSort(T *A, int low, int high, bool parallel);
template <typename T>
void PrintArray(T *A, const int n);
template <typename T>
inline int CompareValues(T a, T b);
template <typename T>
inline bool Equal(T a, T b, double tolerance=0.01);
template <typename T>
string TestOrdering(T *A, size_t n);


int main(int argc, char *argv[])
{
  // 1 million   =  1000000
  // 10 million  =  10000000
  // 100 million =  100000000
  // 1 billion   =  1000000000
  int n; // Size of the array

  if(argc == 1)
    {
      cerr << "\nNo input size given" << endl;
      cerr << "Defaulting to array size 1,000,000\n" << endl;
      n = 1000000;
    }
  else
    n = atoi(argv[1]);

  cout << "\nAllocating Memory & Building Arrays...";
  cout.flush();
  int *ints = new int[n];
  int *ints_seq = new int[n];
  int *ints_std_sort = new int[n];
  //cout << "\tDONE" << endl;
  
  cout << setprecision(4) << fixed;

  // Initiates the seed for rand() based off of current system time
  srand(time(NULL));
  
  // Builds the array
  double makeStart = omp_get_wtime();
  #pragma omp parallel for
  for(int i = 0; i < n; ++i)
    {
      ints_seq[i] = rand() % 113;
      ints[i] = ints_seq[i];
      ints_std_sort[i] = ints[i];
    }
  cout << "\tDONE" << endl;
  cout << "\tTime to build array: \t\t" << omp_get_wtime() - makeStart
       <<  " seconds" << endl;

  cout << "\nBefore sorting...\n";
  cout << "-----------------\n";

  cout << "Testing Array of Integers: \t\t";
  cout.flush();
  cout << TestOrdering(ints, n) << endl;

  // Sequential
  cout << "Performing Sequential Quicksort... ";
  cout.flush();
  double startTime = omp_get_wtime();
  QuickSort(ints_seq, n, sizeof(ints[0]), CompareValues, false);
  cout << "\tDONE" << endl;
  cout << "\tTime to sort: \t\t\t" << omp_get_wtime() - startTime
       << " seconds" << endl;

  // Parallel
  cout << "Performing Parallel Quicksort... ";
  cout.flush();
  startTime = omp_get_wtime();
  QuickSort(ints, n, sizeof(ints[0]), CompareValues);
  cout << "\tDONE" << endl;
  cout << "\tTime to sort: \t\t\t" << omp_get_wtime() - startTime
       << " seconds" << endl;

  // std::sort()
  cout << "Performing std::sort()... ";
  cout.flush();
  startTime = omp_get_wtime();
  sort(ints_std_sort, ints_std_sort + n);
  cout << "\t\tDONE" << endl;
  cout << "\tTime to sort: \t\t\t" << omp_get_wtime() - startTime
       << " seconds" << endl;

  // Making sure it sorted correctly
  cout << "\nAfter sorting...\n";
  cout << "----------------\n";
  
  cout << "Testing Array of Integers: \t\t";
  cout.flush();
  cout << TestOrdering(ints, n) << endl;

  // Freeing memory
  delete [] ints;
  delete [] ints_seq;
  delete [] ints_std_sort;

  return 0;
}

// Generalized QuickSort function that sorts the entire array of data
// A          = Generalized array of data
// n          = Number of elements in the array
// el_size    = Size of an individual element
// Comparison = Pointer to the function to be used to comapre the data
template <typename T>
void QuickSort(T *A, int n, size_t size, int (*Comparison)(T,T), bool parallel)
{
  QSort(A, 0, n-1, parallel);
}

// Templated function of quicksort so that it will work for any type of array
// as long as the comparison function is defined for the data being passed
// A          = Generalized array of data
// low        = Smallest element to compare/sort with
// high       = Highest element to compare/sort with
// Adapted from: http://cs.fit.edu/~pkc/classes/writing/hw15/song.pdf
template <typename T>
void QSort(T *A, int low, int high, bool parallel)
{
  int low_temp = low;
  int high_temp = high;
  T midPoint = A[(low+high)/2];

  while(low_temp <= high_temp)
    {
      while(A[low_temp] < midPoint)
	++low_temp;
      while(A[high_temp] > midPoint)
	--high_temp;
      if(low_temp <= high_temp)
	{
	  swap(A[low_temp], A[high_temp]);
	  ++low_temp;
	  --high_temp;
	}
    }

  #pragma omp parallel sections if(parallel == true)
  {
    #pragma omp section
    if(low < high_temp)
      QSort(A, low, high_temp, parallel);
    #pragma omp section
    if(low_temp < high)
      QSort(A, low_temp, high, parallel);
  }
  return;
}

// Templated Partition function which partitions the passed array. It determins
// the "midpoint" based on the high and low which are passed in, and it is
// based off of the comparision function that is passed into the QSort
// function, which has to be tailored to the data that is going to be sorted
// A          = Generalized array of data
// low        = Smallest element to compare/sort with
// high       = Highest element to compare/sort with
// n          = Number of elements in the array
// Comparison = Pointer to the function to be used to comapre the data
template <typename T>
int Partition(T *A, int low, int high, size_t n, int (*Comparison)(T,T))
{
  T pivot = A[high]; // The "pivot" point
  int i = low - 1;
  int private_i = 0;
  
  #pragma omp parallel for reduction(+:private_i) if(1==0)
  for(int j = low; j < high; ++j)
    {
      // Swaps if the comparision function returns -1, which denotes
      // that it's opposite of the preferred sorting
      if(Comparison(A[j], pivot) == -1)
	{
	  ++private_i;
	  swap(A[private_i + i], A[j]);
	}
    }

  swap(A[i + 1 + private_i], A[high]);
  return i + 1 + private_i;
}

// Prints out each element of the generalized array line by line.
// A = Generalized array of data
// n = Number of elements in the array
template <typename T>
void PrintArray(T *A, const int n)
{
  cout << fixed << setprecision(2);
  for(int i = 0; i < n; ++i)
    cout << A[i] << endl;
    
  return;
}

// Generalized comparision function for int, double, float, long, etc. 
// To be used for nuemerical data types. Defaults to ascending order 
// if not defined in the function call.
// a = General standard numerical data type
// b = General standard numerical data type
template <typename T>
inline int CompareValues(T a, T b)
{
  if(a > b)
    return 1;
  else
    return -1;
}

// Compares two decimal values (float, double, etc) and checks to see if the
// difference between the two is smaller than some preset tolerance. If it is,
// then they are considered equal. The main point of this function is to
// check to see if two decimal numbers are the same without worrying about
// numerical errors
// a         = General standard numerical data type
// b         = General standard numerical data type
// tolerance = Minimum limit of difference to be considered the "same"
template <typename T>
inline bool Equal(T a, T b, double tolerance)
{
  if(abs(a - b) < tolerance)
    return true;
  else
    return false;
}

// Tests the ordering of the given array to see if it was ordered correctly.
// A = Generalized array of data
// n = Number of elements in the array
template <typename T>
string TestOrdering(T *A, size_t n)
{
  for(int i = 0; i < n-1; ++i)
    {
      if(A[i] > A[i+1])
	return "FAILED";
    }
  return "PASSED";
}
