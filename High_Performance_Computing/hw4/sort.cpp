#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <mpi.h>
#include "bitonic.h"
#include "bitonic_MPI.h"

using namespace std;

int main(int argc, char *argv[])
{
  MPI::Init(argc, argv);
  
  int size = 10;
  double wTime_seq = 0;
  double wTime_par = 0;
  double wTime = 0;
  // Sets the MPI settings that will be used
  MPI_Settings mpi;
  mpi.node_id = MPI::COMM_WORLD.Get_rank();
  mpi.nodes = MPI::COMM_WORLD.Get_size();
  mpi.comm = MPI::COMM_WORLD.Dup();
  mpi.type = MPI::INT; // Type sorted with MPI
  mpi.check = false;   // Checks the ordering after sorting

  if(mpi.nodes % 2 != 0)
    {
      if(mpi.node_id == mpi.root)
	cerr << "Error: Need nodes to be a power of 2!\nExiting..." << endl;
      MPI::Finalize();
      return 0;
    }
  
  if(argc == 1)
    size = 1000;
  else
    size = atoi(argv[1]);

  vector<int> A(size);

  if(mpi.node_id == mpi.root)
    {
      cout << "Running for " << mpi.nodes << " nodes with N = "
	   << size << " on each node" << endl;
    }

  // Does this hurt since you're vegetarian?
  srand(0xDEADBEEF + mpi.node_id);

  for(int i = 0; i < A.size(); ++i)
    A[i] = rand();

  if(mpi.node_id == mpi.root)
    cout << "Performing Sequential sort... ";
  // Performs sequential (local) bitonic sort
  wTime_seq = MPI::Wtime();
  if(size > 1)
    BitonicSort(A);
  wTime_seq = MPI::Wtime() - wTime_seq;
  if(mpi.node_id == mpi.root)
    cout << "DONE" << endl;

  if(mpi.node_id == mpi.root)
    cout << "Performing Parallel sort..... ";
  // Performs the sorting between the nodes
  wTime_par = MPI::Wtime();
  if(mpi.nodes > 1)
    Bitonic_MPI(A, mpi);
  wTime_par = MPI::Wtime() - wTime_par;
  if(mpi.node_id == mpi.root)
    cout << "DONE" << endl;

  wTime_par += wTime_seq;

  // Gets the maximum time that a node took (sequential + parallel)
  mpi.comm.Reduce(&wTime_par, &wTime, 1, MPI::DOUBLE, MPI::MAX, 0);

  if(mpi.node_id == mpi.root)
    cout << "Max Execution Time: " << wTime_par << " seconds" << endl;
  
  // Checks the ordering if defined as true above
  if(mpi.check == true)
    CheckOrdering(A, mpi);

  MPI::Finalize();

  return 0;
}
