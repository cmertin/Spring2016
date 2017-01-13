#include <iostream>
#include <cstdlib>
#include <vector>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{
  MPI::Init(argc, argv);
  int size = 100;
  int node_id = MPI::COMM_WORLD.Get_rank();
  int nodes = MPI::COMM_WORLD.Get_size();
  vector<int> A(size);
  vector<int> B(size);
  int pivot;

  srand(time(NULL) * node_id);
  for(int i = 0; i < size; ++i)
    A[i] = rand();

  if(node_id == 0)
    pivot = A.size() / 2;

  MPI::COMM_WORLD.Bcast(&pivot, 1, MPI::INT, 0);

  // Parititioning
  int count = 0;
  for(int i = 0; i < A.size(); ++i)
    {
      if(A[i] < pivot) // Counts the number less than the pivot
	++count;
    }
  int off_lo = 0;
  int off_hi = count;
  for(int i = 0; i < A.size(); ++i)
    {
      if(A[i] < pivot)
	B[off_lo++] = A[i];
      else
	B[off_hi++] = A[i];
    }

  int newNodes = nodes >> 1;
  int is_lo = node_id < newNodes;
  int partner = (node_id + newNodes) % nodes;
  
  

  MPI::Finalize();
  return 0;
}
