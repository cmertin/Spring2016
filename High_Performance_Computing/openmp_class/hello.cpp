#include <iostream>
#include <cstdlib>
#include <omp.h>

using namespace std;

int main(int argc, char *argv[])
{
  if(argc < 2)
    {
      cerr << "Error. Call as: \"" << argv[0] << " num_threads\"" << endl;
      return -1;
    }
  int numThreads = atoi(argv[1]);

#pragma omp parallel num_threads(numThreads)
  {
    int tid = omp_get_thread_num();
    cout << "Hello world from thread " << tid << endl;
  }

  return 0;
}
