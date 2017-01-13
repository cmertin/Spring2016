#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <omp.h>
#include <algorithm>
#include <mpi.h>
#include "bsearch.h"

using namespace std;

int main(int argc, char *argv[])
{
  MPI::Init(argc, argv);
  int size;
  if(argc == 1)
      size = 1000000;
  else
    size = atoi(argv[1]);

  cout << setprecision(5) << fixed;
  parallelSettings<int> parallel;

  int numKeys = 4;
  parallel.keys.resize(numKeys);
  parallel.node_id = MPI::COMM_WORLD.Get_rank();
  parallel.nodes = MPI::COMM_WORLD.Get_size();
  parallel.numThreads = omp_get_max_threads();
  parallel.comm = MPI::COMM_WORLD.Dup();
  parallel.numPerNode = size/parallel.nodes;
  double startTime = 0;
  double endTime = 0;
  double maxTime = 0;
  double minTime = 0;
  vector<int> A(parallel.numPerNode);
  vector<int>::iterator itrBegin = A.begin();
  vector<int>::iterator itrEnd = A.end();

  if(parallel.node_id == 0)
    {
      cout << "Running for N = " << size << endl;
      cout << "Number of nodes = " << parallel.nodes << endl;
    }
    
  srand(time(NULL) * parallel.node_id);

  // Get the keys on the first node
  // They will be broadcasted in the function
  if(parallel.node_id == 0)
    {
      for(int i = 0; i < parallel.keys.size(); ++i)
	parallel.keys[i] = rand();
    }

  
  if(parallel.node_id == 0)
    {
      cout << "Building array...... ";
      cout.flush();
    }
  
  // Builds the values in the arrays
  #pragma parallel for num_threads(parallel.numThreads)
  for(auto itr = itrBegin; itr != itrEnd; ++itr)
    *itr = rand();

  if(parallel.node_id == 0)
    cout << "DONE" << endl;

  // Sorts the arrays using std::sort
  if(parallel.node_id == 0)
    {
      cout << "Sorting array....... ";
      cout.flush();
    }
  startTime = MPI::Wtime();
  sort(itrBegin, itrEnd);
  endTime = MPI::Wtime();
  startTime = endTime - startTime;
  parallel.comm.Reduce(&startTime, &minTime, 1, MPI::DOUBLE, MPI::MIN, 0);
  parallel.comm.Reduce(&startTime, &maxTime, 1, MPI::DOUBLE, MPI::MAX, 0);
  if(parallel.node_id == 0)
    cout << "DONE (" << minTime << " s, " << maxTime << " s)" << endl;

  // Performs the binary search
  if(parallel.node_id == 0)
    {
      cout << "Searching........... ";
      cout.flush();
    }
  startTime = MPI::Wtime();
  vector<int> positions = BinarySearch(itrBegin, itrEnd, parallel);
  endTime = MPI::Wtime();
  startTime = endTime - startTime;
  parallel.comm.Reduce(&startTime, &minTime, 1, MPI::DOUBLE, MPI::MIN, 0);
  parallel.comm.Reduce(&startTime, &maxTime, 1, MPI::DOUBLE, MPI::MAX, 0);
  if(parallel.node_id == 0)
    cout << "DONE (" << minTime << " s, " << maxTime << " s)" << endl;

  MPI::Finalize();
  return 0;
}


