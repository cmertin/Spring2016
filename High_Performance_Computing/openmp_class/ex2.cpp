#include <iostream>
#include <cstdlib>
#include <omp.h>

using namespace std;

double fn(double x);

int main(int argc, char *argv[])
{
  if(argc < 2)
    {
      cerr << "Usage: " << argv[0] << " #" << endl;
      return -1;
    }

  int n = atoi(argv[1]);
  int numThreads = omp_get_max_threads();

  double x, pi, step, t1, t2, sum=0.0;
  double *sumP = new double[numThreads];

  for(int i = 0; i < numThreads; ++i)
    sumP[i] = 0.0;

  step = 1.0/n;

  t1 = omp_get_wtime();
  for(int i = 0; i < n; ++i)
    {
      x += step;
      sum += fn(x);
    }
  pi = step * sum;
  t2 = omp_get_wtime();

  cout << "Sequential Pi: " << pi << endl;
  cout << "Sequential Time: " << t2-t1 << endl << endl;
  
  t1 = omp_get_wtime();
  #pragma omp parallel for
  for(int i = 0; i < n; ++i)
    {
      double xP = (i + 0.5) * step;
      int tid = omp_get_thread_num();
      sumP[tid] += fn(xP);
    }
  pi = step * sum;
  t2 = omp_get_wtime();

  cout << "Parallel Pi: " << pi << endl;
  cout << "Parallel Time: " << t2 - t1 << endl << endl;

  sum = 0.0;

  // Local sum
  t1 = omp_get_wtime();
#pragma omp parallel for reduction(+:sum)
  for(int i = 0; i < n; ++i)
    {
      double xL = (i + 0.5) * step;
      sum += fn(xL);
    }
  pi = step * sum;
  t2 = omp_get_wtime();

// #pragma omp atomic
// pi += sum * step;

cout << "Reduction Pi: " << pi << endl;
cout << "Reduction Time: " << t2 - t1 << endl;


  return 0;
}

inline double fn(double x)
{
  return 4 / (1.0 + x*x);
}
