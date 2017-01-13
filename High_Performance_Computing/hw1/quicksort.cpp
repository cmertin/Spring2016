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

using namespace std;

struct Point
{
  double x;
  double y;
};

template <typename T>
void QuickSort(T *A, size_t num, size_t el_size, int (*Comparison)(T,T));
template <typename T>
void QSort(T *A, int low, int high, size_t n, int (*Comparison)(T,T));
template <typename T>
int Partition(T *A, int low, int high, size_t n, int (*Comparison)(T,T));
void PrintPoints(Point *A, const int n);
template <typename T>
void PrintArray(T *A, const int n);
int ComparePoints(Point a, Point b);
template <typename T>
int CompareValues(T a, T b);
template <typename T>
bool Equal(T a, T b, double tolerance=0.01);
template <typename T>
string TestOrdering(T *A, size_t n);
string TestPointsOrdering(Point *A, size_t n);


int main()
{
  string ordering = "ascending";
  int n = 1000; // Size of the array
  cout << "Array Size: ";
  cin >> n;
  int *ints = new int[n];
  long *longs = new long[n];
  float *floats = new float[n];
  double *doubles = new double[n];
  Point *points = new Point[n]; // Array of points

  // Initiates the seed for rand() based off of current system time
  srand(time(NULL));

  // Builds the array of points with random x and y values
  for(int i = 0; i < n; ++i)
    {
      points[i].x = (rand() % 113) + (rand() % 47) * M_PI;
      points[i].y = (rand() % 101) + (rand() % 59) * M_PI;

      ints[i] = points[i].y;
      longs[i] = points[i].y;
      floats[i] = points[i].y;
      doubles[i] = points[i].y;
    }

  cout << "Before sorting...\n";
  cout << "-----------------\n";

  cout << "Testing Array of Integers: \t"
       << TestOrdering(ints, n) << endl;
  cout << "Testing Array of Longs: \t"
       << TestOrdering(longs, n) << endl;
  cout << "Testing Array of Floats: \t"
       << TestOrdering(floats, n) << endl;
  cout << "Testing Array of Doubles: \t"
       << TestOrdering(doubles, n) << endl;
  cout << "Testing Array of Points: \t"
       << TestPointsOrdering(points, n) << endl;

  QuickSort(ints, n, sizeof(ints[0]), CompareValues);
  QuickSort(longs, n, sizeof(longs[0]), CompareValues);
  QuickSort(floats, n, sizeof(floats[0]), CompareValues);
  QuickSort(doubles, n, sizeof(doubles[0]), CompareValues);
  QuickSort(points, n, sizeof(points[0]), ComparePoints);

  cout << "\nAfter sorting...\n";
  cout << "----------------\n";
  
  cout << "Testing Array of Integers: \t"
       << TestOrdering(ints, n) << endl;
  cout << "Testing Array of Longs: \t"
       << TestOrdering(longs, n) << endl;
  cout << "Testing Array of Floats: \t"
       << TestOrdering(floats, n) << endl;
  cout << "Testing Array of Doubles: \t"
       << TestOrdering(doubles, n) << endl;
  cout << "Testing Array of Points: \t"
       << TestPointsOrdering(points, n) << endl;
  
  return 0;
}

// Generalized QuickSort function that sorts the entire array of data
// A          = Generalized array of data
// n          = Number of elements in the array
// el_size    = Size of an individual element
// Comparison = Pointer to the function to be used to comapre the data
template <typename T>
void QuickSort(T *A, size_t n, size_t el_size, int (*Comparison)(T,T))
{
  QSort(A, 0, n-1, n, Comparison);
}

// Templated function of quicksort so that it will work for any type of array
// as long as the comparison function is defined for the data being passed
// A          = Generalized array of data
// low        = Smallest element to compare/sort with
// high       = Highest element to compare/sort with
// n          = Number of elements in the array
// Comparison = Pointer to the function to be used to comapre the data
template <typename T>
void QSort(T *A, int low, int high, size_t n, int (*Comparison)(T,T))
{
  if(low < high)
    {
      // Defines the midpoint of the current segment
      int mid = Partition(A, low, high, n, Comparison);
      QSort(A, low, mid-1, n, Comparison);
      QSort(A, mid+1, high, n, Comparison);
    }
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
  for(int j = low; j < high; ++j)
    {
      // Swaps if the comparision function returns -1, which denotes
      // that it's opposite of the preferred sorting
      if(Comparison(A[j], pivot) == -1)
	{
	  ++i;
	  swap(A[i],A[j]);
	}
    }
  swap(A[i+1],A[high]);
  return i + 1;
}

// Prints the array of points line by line with fixed 2 decimals. Adds extra
// spaces to make it look "nice" and to verify results.
// A = Generalized array of data
// n = Number of elements in the array
void PrintPoints(Point *A, const int n)
{
  cout << fixed << setprecision(2);
  for(int i = 0; i < n; ++i)
    {
      if(A[i].x < 100)
	cout << "( " << A[i].x << ", ";
      else
	cout << '(' << A[i].x << ", ";

      if(A[i].y < 100)
	cout << ' ' << A[i].y << ")\n";
      else
	cout << A[i].y << ")\n";
    }
  return;
}

// Prints out each element of the generalized array line by line.
// A = Generalized array of data
// n = Number of elements in the array
template <typename T>
void PrintArray(T *A, const int n)
{
  cout << fixed << setprecision(2);
  for(int i = 0; i < n; ++i)
    {
      if(A[i] < 100)
	cout << ' ' << A[i] << '\n';
      else
	cout << A[i] << '\n';
    }
  return;
}

// Compares the points and returns 1 if it matches the sort order and -1 if not
// the default sorting order is "ascending" though "descending" can also be
// passed in as well, and it will sort it in descending order instead. It sorts
// it first by the y-values, and if the y-values are the same, then it compares
// based off of the x-values.
// a = Point element to be compared
// b = Point element to be compared
int ComparePoints(Point a, Point b)
{
  if(a.y > b.y && !Equal(a.y, b.y))
    return 1;
  else
    {
      if(Equal(a.y, b.y) && a.x > b.x)
	return 1;
      else
	return -1;
    }
}

// Generalized comparision function for int, double, float, long, etc. 
// To be used for nuemerical data types. Defaults to ascending order 
// if not defined in the function call.
// a = General standard numerical data type
// b = General standard numerical data type
template <typename T>
int CompareValues(T a, T b)
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
bool Equal(T a, T b, double tolerance)
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

// Tests the ordering of a given array of type Points to
// see if they're ordered correctly
// A = Array of type Points
// n = Number of elements in the array
string TestPointsOrdering(Point *A, size_t n)
{
  for(int i = 0; i < n-1; ++i)
    {
      if(A[i].y > A[i+1].y && !Equal(A[i].y, A[i+1].y))
	{
	  return "FAILED";
	}
      else if(Equal(A[i].y, A[i+1].y) && A[i].x > A[i+1].x)
	{
	  return "FAILED";
	}
    }
  return "PASSED";
}
