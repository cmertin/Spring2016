#ifndef BSEARCH_H
#define BSEARCH_H
#include <vector>
#include <omp.h>
#include <mpi.h>

template<typename T>
struct parallelSettings
{
  std::vector<T> keys;      // Keys to search for
  MPI::Intracomm comm;      // MPI Communicator
  int numThreads;           // Number of threads
  int node_id;              // ID of each node
  int nodes;                // Number of nodes
  int numPerNode;           // Number of elements per node
};

template <typename Iterator, typename T>
std::vector<int> BinarySearch(Iterator begin, Iterator end, parallelSettings<T> &par)
{
  // Broadcast the keys to the other nodes
  par.comm.Bcast(&par.keys.front(), par.keys.size(), MPI::INT, 0);

  int numElements = distance(begin, end); // Gets the number of elements
  int numElementsPer = distance(begin,end)/par.numThreads; // Number of elements per OpenMP thread
  bool isFound = false; // If the key is found
  int index = 0; // Indexing for the positions
  std::vector<int> positions(par.keys.size(), -1); // Holds the positions
  std::vector<int> allPositions(par.keys.size() * par.nodes); // Holds all positions

  for(Iterator keys = par.keys.begin(); keys != par.keys.end(); ++keys)
    {
      T key = *keys;
      isFound = false;
      #pragma omp parallel num_threads(par.numThreads) shared(isFound)
      {
	int threadNum = omp_get_thread_num(); // Number of threads
	Iterator beginT = begin + threadNum * numElementsPer; // Beginning element for the thread
	Iterator endT = beginT + numElementsPer; // End element for the thread
	for(Iterator itr = beginT; itr < endT; ++itr)
	  {
	    if(isFound == true)
	      continue;
	    else
	      {
		Iterator middle = beginT + (distance(beginT,endT)/2);
		if(*middle == key)
		  {
		    positions[index] = distance(begin, middle); // Save the found position
		    isFound = true;
                    #pragma omp flush (isFound) // updates isFound for all threads and registers
		  }
		  else if(*middle > key)
		    endT = middle;
		  else
		    beginT = middle + 1;
		  }
	      }
	  }
      ++index;
      }

    // Appends the minimum value for each node since the array is separated
    if(par.node_id > 0)
      {
	for(int i = 0; i < positions.size(); ++i)
	  {
	    if(positions[i] != -1)
	      positions[i] += par.numPerNode * par.node_id;
	  }
      }
    // Takes all of the values from each node to create a vector
    // and gives it to each node
    par.comm.Allgather(&positions.front(), positions.size(), MPI::INT,
		    &allPositions.front(), positions.size(), MPI::INT);

  return allPositions;
}


#endif
