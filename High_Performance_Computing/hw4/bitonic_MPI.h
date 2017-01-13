#ifndef BITONIC_MPI_H
#define BITONIC_MPI_H
#include <iostream>
#include <mpi.h>
#include <vector>
#include "bitonic.h"

struct MPI_Settings
{
  MPI::Intracomm comm; // The MPI Communicator
  MPI::Datatype type;  // The data type to be sorted
  int node_id;         // ID/rank of each node
  int nodes;           // Total number of nodes
  bool check = false;  // Used to check the ordering. Defaults to false.
  int root = 0;        // Defines the "root node" (defaults to 0)
};

// Checks the ordering at the end of the sorting. It gathers all the values
// (in order) to a single node and then checks to make sure the ordering
// is sequential.
template <typename T>
void CheckOrdering(std::vector<T> &A, MPI_Settings &mpi)
{
  std::vector<T> allVals(A.size() * mpi.nodes);
  mpi.comm.Gather(&A[0], A.size(), mpi.type, &allVals[0], A.size(),
		  MPI::INT, mpi.root);
  if(mpi.node_id == mpi.root)
    {
      std::cout << "Checking Ordering............ "; 
      for(int i = 0; i < allVals.size()-1; ++i)
	{
	  if(allVals[i] <= allVals[i + 1])
	    continue;
	  else
	    {
	      std::cout << "FAILED" << std::endl;
	      return;
	    }
	}
      std::cout << "PASSED" << std::endl; 
    }
  return;
}

// Compares against the low values. It sends its pivot/mid value to the other
// node, and then iterates through comparing the received value to the lower
// values in the array to find out if anything needs to be swapped
template <typename T>
void CompareLow(std::vector<T> &A, int otherNode, MPI_Settings &mpi)
{
  int low = 0;
  int high = A.size();
  int pivot = A.size()/2;
  T otherVal = 0;
  T tempVal = A[pivot];
  mpi.comm.Send(&tempVal, 1, mpi.type, otherNode, 0);
  mpi.comm.Recv(&otherVal, 1, mpi.type, otherNode, 1);
  
  while(pivot != low && A[pivot] != otherVal)
    {
      if(A[pivot] < otherVal)
	low = pivot;
      else
	high = pivot;
      pivot = low + (high - low)/2;
      tempVal = A[pivot];

      mpi.comm.Send(&tempVal, 1, mpi.type, otherNode, 0);
      mpi.comm.Recv(&otherVal, 1, mpi.type, otherNode, 1);
    }
 
  if(A[pivot] <= otherVal)
    {
      int size = pivot + 1;
      std::vector<T> temp(size);
      mpi.comm.Send(&A[0], size, mpi.type, otherNode, 0);
      mpi.comm.Recv(&temp[0], size, mpi.type, otherNode, 1);

      for(int i = 0; i < temp.size(); ++i)
	swap(A[i], temp[i]);
    }
  else if(otherVal < A[pivot])
    {
      int size = pivot;
      std::vector<T> temp(size);
      mpi.comm.Send(&A[0], size, mpi.type, otherNode, 0);
      mpi.comm.Recv(&temp[0], size, mpi.type, otherNode, 1);
      
      for(int i = 0; i < temp.size(); ++i)
	swap(A[i], temp[i]);
    }
  return;
}

// Performs the opposite of the CompareLow function. This one compares the
// upper values to the received value
template <typename T>
void CompareHigh(std::vector<T> &A, int otherNode, MPI_Settings &mpi)
{
  int low = -1;
  int high = A.size() - 1;
  int pivot = low + (high - low + 1)/2;
  T otherVal = 0;
  T tempVal = A[pivot];
  mpi.comm.Recv(&otherVal, 1, mpi.type, otherNode, 0);
  mpi.comm.Send(&tempVal, 1, mpi.type, otherNode, 1);
    
  while(pivot != low && pivot != high && A[pivot] != otherVal)
    {
      if(A[pivot] < otherVal)
	low = pivot;
      else
	high = pivot;
      pivot = low + (high - low + 1)/2;
      tempVal = A[pivot];

      mpi.comm.Recv(&otherVal, 1, mpi.type, otherNode, 0);
      mpi.comm.Send(&tempVal, 1, mpi.type, otherNode, 1);
    }
  
  if(otherVal <= A[pivot])
    {
      int size = A.size() - pivot;
      std::vector<T> temp(size);
      mpi.comm.Recv(&temp[0], size, mpi.type, otherNode, 0);
      mpi.comm.Send(&A[pivot], size, mpi.type, otherNode, 1);

      for(int i = 0; i < temp.size(); ++i)
	swap(A[pivot + i], temp[i]);
    }
  else if(A[pivot] < otherVal)
    {
      int size = A.size() - pivot - 1;
      std::vector<T> temp(size);
      mpi.comm.Recv(&temp[0], size, mpi.type, otherNode, 0);
      mpi.comm.Send(&A[pivot + 1], size, mpi.type, otherNode, 1);

      for(int i = 0; i < temp.size(); ++i)
	swap(A[pivot + i], temp[i]);
    }
  return;
}

// This calculates the nodes to send and receive from so that it is done in
// a "bitonic" style. The pseudocode where this code was derived from was
// taken from [http://web.mst.edu/~ercal/387/P3/pr-proj-3.pdf] and adapted
// into C++ code.
template <typename T>
void Bitonic_MPI(std::vector<T> &A, MPI_Settings &mpi)
{
  for(int i = 1; i <= log2(mpi.nodes); ++i)
    {
      int node = (mpi.node_id >> i) & 1;
      for(int j = i-1; j >= 0; --j)
	{
	  int bit1 = (mpi.node_id >> j) & 1;
	  int sendID = mpi.node_id ^ (1 << j);
	  if(node != bit1)
	    CompareLow(A, sendID, mpi);
	  else
	    CompareHigh(A, sendID, mpi);

	  if(A.size() > 1)
	    BitonicSort(A);
	}
    }
  return;
}

#endif
