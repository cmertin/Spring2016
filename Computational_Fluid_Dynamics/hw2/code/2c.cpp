#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <ctime>

using namespace std;

template <typename T>
struct Matrix
{
  vector<T> lower; 
  vector<T> diag;  
  vector<T> upper; 
  int n;           // dimensions of matrix
};

template <typename T>
void TridiagSolve(Matrix<T> A, vector<T> b, vector<T> &x);
void WriteFile(vector<double> x, vector<double> u, string time);

int main()
{
  double h = 0.04;
  double dy = 0.0005;
  double dy2 = dy*dy;
  double dt = 0.0005;//0.00333;
  int n = h/dy;
  double nu = 0.000217;
  double rho = 800.0;
  double dpdx = 2000.0;
  double beta = 1/rho * dpdx;
  int num_itr = 1/dt;
  Matrix<double> A;
  A.diag.resize(n);
  A.upper.resize(n);
  A.lower.resize(n);
  A.n = n;
  vector<double> b(n);
  vector<double> u(n);
  vector<double> x(n);
  double time = 0.0;
  int count = 0;
  int write_count = 0.2/dt;
  chrono::time_point<chrono::system_clock> start,end;
  cout << "Time stepping: dt = " << dt << endl;
  start = chrono::system_clock::now();
  // Fills the matrix
  for(int i = 0; i < A.n; ++i)
    {
      double alpha = nu/(2 * dy2);
      A.diag[i] = (1.0/dt + 2*alpha);
      if(i < n-1)
	A.upper[i] = -alpha;
      if(i > 0)
	A.lower[i] = -alpha;

    }
  x[0] = 0.0;
  A.diag[0] = 1; // Set for B.C.
  A.diag[A.n-1] = 1; // Set for B.C.
  A.upper[0] = 0; // Set for B.C.
  A.lower[A.n-1] = 0; // Set for B.C.

  u[0] = 40.0;
  for(int i = 1; i < u.size(); ++i)
    {
      u[i] = 0.0;
      x[i] = x[i-1] + dy;
    }

  WriteFile(x, u, to_string(time));
  for(int i = 0; i <= num_itr; ++i)
    {
      time += dt;
      double alpha = nu/(2 * dy2);
      // Build b
      for(int i = 0; i < b.size(); ++i)
	b[i] = (1.0/dt - 2*alpha)*u[i] + alpha*u[i-1] + alpha*u[i+1] - beta;
      // Impose B.C. 
      b[0] = 40.0;
      b[A.n-1] = 0.0;
      TridiagSolve(A, b, u);
      if(count == write_count-1)
	{
	  WriteFile(x, u, to_string(time));
	  count = 0;
	}
      else
	++count;
    }
  end = chrono::system_clock::now();
  chrono::duration<double> runtime = end-start;

  cout << "Run time: " << runtime.count() << " seconds" << endl;
  
  return 0;
}


// Uses the Thomas Algorithm
// Assumes A is positive semi-definite and tridiagonal
template <typename T>
void TridiagSolve(Matrix<T> A, vector<T> b, vector<T> &x)
{
  // Forward Elimination
  A.upper[0] = A.upper[0]/A.diag[0];
  b[0] = b[0]/A.diag[0];
  for(int i = 1; i < A.n; ++i)
    {
      if(i < A.n-1)
	A.upper[i] = A.upper[i]/(A.diag[i] - A.lower[i]*A.upper[i-1]);
      b[i] = (b[i] - A.lower[i]*b[i-1])/(A.diag[i] - A.lower[i]*A.upper[i-1]);
    }

  // Backward Substitution
  x[A.n-1] = b[A.n-1];
  for(int i = A.n-2; i >= 0; --i)
    x[i] = b[i] - A.upper[i]*x[i+1];

  return;
}

void WriteFile(vector<double> x, vector<double> u, string time)
{
  stringstream ss;
  ss << time[0] << time[1] << time[2] << time[3] << time[4];// << time[5] << time[6];
  string true_time = ss.str();
  cout << true_time << " seconds" << endl;
  string outputfile = true_time + "_sec.dat";
  string sys_cmd = "./plot_2.py " + true_time;
  ofstream outfile;
  outfile.open(outputfile.c_str());
  for(int i = 0; i < u.size(); ++i)
    {
      outfile << x[i] << ',' << u[i] << endl;
    }
  outfile.close();

  system(sys_cmd.c_str());

  return;
}
