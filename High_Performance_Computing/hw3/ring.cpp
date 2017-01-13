#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
  MPI::Init(argc, argv);
  int loops;
  int total = 0;
  double wTime;
  int node_id = MPI::COMM_WORLD.Get_rank();
  int nodes = MPI::COMM_WORLD.Get_size();
  int bigN = 10000000;
  vector<int> largeVector(bigN);
  int size = sizeof(largeVector);

  // Default to 1 loop
  if(argc == 1)
    {
      loops = 1;
      if(node_id == 0)
	cerr << "Defaulting to a single loop." << endl;
    }
  
  // Number of loops around the ring
  else if(argc == 2)
    loops = atoi(argv[1]);

  // Close program for incorrect arguments
  else
    {
      cerr << "Only allows for zero or one arguments. Exiting..." << endl;
      return -1;
    }

  // Initialize time
  if(node_id == 0)
    {
      cout << "Number of nodes: " << nodes << endl;
      cout << "Number of loops: " << loops << endl;
      wTime = MPI::Wtime();
    }

  // Send around the ring
  for(int i = 0; i < loops; ++i)
    {
      for(int j = 1; j < nodes; ++j)
	{
	  total += node_id;
	  if(node_id == j-1)
	    MPI::COMM_WORLD.Send(&total, 1, MPI::INT, j, 0);
	  if(node_id == j)
	    MPI::COMM_WORLD.Recv(&total, 1, MPI::INT, j-1, 0);
	}
      total += node_id;
      if(node_id == nodes-1)
	MPI::COMM_WORLD.Send(&total, 1, MPI::INT, 0, 0);
      if(node_id == 0)
	MPI::COMM_WORLD.Recv(&total, 1, MPI::INT, nodes-1, 0);
    }
  
  // Finialize time
  if(node_id == 0)
    {
      wTime = MPI::Wtime() - wTime;
      cout << "Total: " << total << endl;
      cout << "Run time: " << wTime << " seconds" << endl;
      cout << "Average time: " << wTime/(nodes*loops) << " seconds" << endl;
      cout << "Number of Sends: " << nodes * loops << endl;

      cout << "\n\nTesting for large array" << endl;
      cout <<     "-----------------------" << endl;
      cout << "N = " << bigN << endl;
    }

  wTime = MPI::Wtime();
  for(int i = 1; i < nodes; ++i)
    {
      if(node_id == i-1)
	MPI::COMM_WORLD.Send(&largeVector.front(), bigN, MPI::INT, i, 0);
      if(node_id == i)
	MPI::COMM_WORLD.Recv(&largeVector.front(), bigN, MPI::INT, i-1, 0);
    }
  wTime = MPI::Wtime() - wTime;

  if(node_id == 0)
    {
      cout << "Average time to send the large array: "
	   << wTime / nodes << " seconds" << endl;
      cout << "Size of the array: " << size * bigN << " bytes" << endl;
      cout << "bytes/second: " << (size * bigN)/(wTime/nodes) << endl;
    }


  

  MPI::Finalize();
  return 0;
}
