### BitonicSort
Files: `bitonic.h bitonic_MPI.h sort.cpp`
This implements the paralellel and serial versions of bitonic sort. `bitonic.h`
contains the proper code and function calls for sorting the data sequentially,
and `bitonic_MPI.h` is for sorting the data between the nodes. They were made
separate so that only the serial version can be used if desired. `bitonic_MPI.h`
requires `bitonic.h` to run as it calls functions from that file as well.

To run the code, first run the `make` command to make the file and then run it
by `mpirun -np nodes ./sort N` where `nodes` is the number of nodes that you're
going to use, and `N` is the size of the array on each node. Note: To run the
code, the number of nodes must be a power of 2 or else it will fail. I haven't
figured out how to get it to work on any arbitrary number of nodes as of yet.

#### Bugs/Features:
* Doesn't work for number of Nodes not a power of 2