#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char *argv[])
{
  MPI::Init(argc, argv);
  // MPI::Init(&argc, &argv);
  // MPI_Init(&argc, &argv);
  // MPI_Comm comm = MPI_COMM_WORLD;
  int node_id = MPI::COMM_WORLD.Get_rank();
  int nodes = MPI::COMM_WORLD.Get_size();

  cout << "I am from processor " << node_id << " of " << nodes << endl;
  MPI::COMM_WORLD.Barrier();

  int color = nodes % 2;
  // MPI_Comm_split(MPI_COMM_WORLD, color, node_id, &comm_new);

  cout << "Round 2: I am from processor " << node_id << " of " << nodes << endl;

  MPI::Finalize();
  return 0;
}
