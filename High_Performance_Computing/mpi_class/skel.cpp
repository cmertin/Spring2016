#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{
  MPI::Init(argc, argv);
  int node_id = MPI::COMM_WORLD.Get_rank();
  int nodes = MPI::COMM_WORLD.Get_size();

  MPI::Finalize();
  return 0;
}
