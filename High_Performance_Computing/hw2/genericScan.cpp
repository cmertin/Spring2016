#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <omp.h>
#include "vector.h"

using namespace std;

template <typename T>
void GenericScan(vector<T> &A, T (*Oper)(T &a1, T &a2), bool parallel=true);
template <typename T>
inline T Add(T &a1, T &a2);
template <typename T>
inline T Subtract(T &a1, T &a2);
template <typename T>
inline T Multiply(T &a1, T &a2);
template <typename T>
inline T Divide(T &a1, T &a2);
template <typename T>
void PrintArray(vector<T> &A);

int main(int argc, char *argv[])
{
  int n;

  if(argc == 1)
    {
      cerr << "\nDefaulting to an array of size 1,000,000\n" << endl;
      n = 1000000;
    }
  else
    n = atoi(argv[1]);

  srand(time(NULL));

  if(false)
    {
      cout << setprecision(5) << fixed;
      cout << "Testing cases for 1D and 3D doubles (300,000,000)" << endl;
      int threeM = 300000000;
      cout << "Allocating Array......................... ";
      cout.flush();

      vector<double> oneD(threeM);
      vector<double> oneD_seq(threeM);
      vector<Vector<double> > threeD(threeM);
      vector<Vector<double> > threeD_seq(threeM);
      double startTime = 0;

      cout << "DONE" << endl;

      cout << "Building Array........................... ";
      cout.flush();

      #pragma parallel for
      for(int i = 0; i < threeM; ++i)
	{
	  oneD[i] = (rand() % 997) * M_PI;
	  oneD_seq[i] = oneD[i];
	  threeD[i].SetX((rand() % 911) * M_PI);
	  threeD[i].SetY((rand() % 491) * M_PI);
	  threeD[i].SetZ((rand() % 757) * M_PI);

	  threeD_seq[i].SetX(threeD[i].GetX());
	  threeD_seq[i].SetY(threeD[i].GetX());
	  threeD_seq[i].SetZ(threeD[i].GetX());
	}
      cout << "DONE" << endl;

      cout << "Performing 1D Sequential Scan............ ";
      cout.flush();
      startTime = omp_get_wtime();
      GenericScan(oneD_seq, Add, false);
      cout << "DONE (" << omp_get_wtime() - startTime << " seconds)" << endl;

      cout << "Performing 1D Parallel Scan.............. ";
      cout.flush();
      startTime = omp_get_wtime();
      GenericScan(oneD, Add);
      cout << "DONE (" << omp_get_wtime() - startTime << " seconds)" << endl;

      cout << "Performing 3D Sequential Scan............ ";
      cout.flush();
      startTime = omp_get_wtime();
      GenericScan(threeD_seq, Add, false);
      cout << "DONE (" << omp_get_wtime() - startTime << " seconds)" << endl;

      cout << "Performing 3D Parallel Scan.............. ";
      cout.flush();
      startTime = omp_get_wtime();
      GenericScan(threeD, Add);
      cout << "DONE (" << omp_get_wtime() - startTime << " seconds)" << endl;

    }
  else
    {
      cout << "Allocating Array.......................... ";
      cout.flush();
      vector<Vector<double> > A(n);
      vector<Vector<double> > A_seq(n);

      cout << "DONE" << endl;

      cout << "Building the Array........................ ";
      cout.flush();
      #pragma parallel for
      for(int i = 0; i < n; ++i)
	{
	  A_seq[i].SetX((rand() % 911) * M_PI);
	  A_seq[i].SetY((rand() % 491) * M_PI);
	  A_seq[i].SetZ((rand() % 757) * M_PI);

	  A[i].SetX(A_seq[i].GetX());
	  A[i].SetY(A_seq[i].GetY());
	  A[i].SetZ(A_seq[i].GetZ());
	}

      cout << "DONE" << endl;
      cout << setprecision(5) << fixed;

      // Sequential Scan
      double startTime_seq = omp_get_wtime();
      cout << "Performing Sequential Scan Operation...... ";
      cout.flush();
      GenericScan(A_seq, Add, false);
      cout << "DONE (" << omp_get_wtime() - startTime_seq
	   << " seconds)" << endl;

      // Parallel Scan
      double startTime = omp_get_wtime();
      cout << "Performing Parallel Scan Operation........ ";
      cout.flush();
      GenericScan(A, Add);
      cout << "DONE (" << omp_get_wtime() - startTime << " seconds)" << endl;
    }
  return 0;
}

template <typename T>
void GenericScan(vector<T> &A, T (*Oper)(T &a1, T &a2), bool parallel)
{
  // Upsweep
  for(int i = 0; i < log2(A.size())-1; ++i)
    {
      #pragma omp parallel for if(parallel && A.size() > 10e6)
      for(int j = 0; j < A.size(); j+=(1<<(i+1)))
	{
	  int m = j + (1<<(i+1)) - 1;
	  if(m < A.size())
	    A[j+(1<<(i+1))-1] = Oper(A[j+(1<<i)-1], A[j+(1<<(i+1))-1]);
	}
    }

  // Downsweep
  for(int i = log2(A.size())-1; i > 0; --i)
    {
      #pragma omp parallel for if(parallel && A.size() > 10e6)
      for(int j = 0; j < A.size(); j+=(1<<i))
	{
	  int m = j + (1<<i) + (i-1);
	  if(m < A.size())
	    A[j+(1<<i)+(i-1)] = Oper(A[j+(1<<i)-1], A[j+(1<<i)+(i-1)]);
	}
    }

  return;
}

template <typename T>
inline T Add(T &a1, T &a2)
{
  return a1 + a2;
}

template <typename T>
inline T Subtract(T &a1, T &a2)
{
  return a1 - a2;
}

template <typename T>
inline T Multiply(T &a1, T &a2)
{
  return a1 * a2;
}

template <typename T>
inline T Divide(T &a1, T &a2)
{
  return a1 / a2;
}

template <typename T>
void PrintArray(vector<T> &A)
{
  for(int i = 0; i < A.size(); ++i)
    cout << A[i] << endl;
  return;
}
