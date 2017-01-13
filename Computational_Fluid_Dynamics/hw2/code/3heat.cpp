#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

template <typename T>
struct Matrix
{
  vector<T> lower; // Lower diagonal
  vector<T> diag;  // Diagonal
  vector<T> upper; // Upper diagonal
  int n;           // Dimensions of matrix
};

template <typename T>
struct Dimensions
{
  vector<T> x;
  vector<T> y;
};

template <typename T>
void TridiagSolve(Matrix<T> A, vector<T> b, vector<T> &x);
template <typename T>
void WriteFile(Dimensions<T> pos, vector<double> u, string time);
bool WriteTime(double time, vector<double> writeTimes, double tol=0.01);
template <typename T>
bool SteadyState(Dimensions<T> pos, vector<double> T1, vector<double> T2, double steadyState);

int main()
{
  double Lx = 0.3;
  double Ly = 0.4;
  int Nx = 61;
  int Ny = 81;
  double dx = Lx/(Nx-1);
  double dy = Ly/(Ny-1);
  double diff = 1.1234 * pow(10,-4); // Diffusivity
  double dt = 0.1;
  double dt_2 = dt/2;
  int n = Nx * Ny;
  int num_itr = 40/dt;
  // I.C. and B.C.
  double T0 = 0.0;
  double T1 = 40.0;  // Bottom of the bar
  double T2 = 0.0;   // Left side of bar
  double T3 = 10.0;  // Top of bar
  double T4 = 0.0;   // Right side of bar
  double steadyState = 0.0001;
  Dimensions<double> pos;
  Matrix<double> Ax;
  Matrix<double> Ay;
  Ax.diag.resize(n);
  Ay.diag.resize(n);
  Ax.upper.resize(n);
  Ay.upper.resize(n);
  Ax.lower.resize(n);
  Ay.lower.resize(n);
  Ax.n = n;
  Ay.n = n;
  vector<double> b(n);
  vector<double> temp(n);
  vector<double> old_temp(n);
  vector<double> x(n);
  double time = 0.0;
  vector<double> writeTimes = {0.5, 10.0, 20.0, 40.0};

  /*
  writeTimes.push_back(0.0);
  for(int i = 1; i < num_itr; ++i)
    writeTimes.push_back(writeTimes[i-1] + dt);
  */
  // Fills for I.C.'s and B.C.'s
  for(int i = 0; i < Ny; ++i)
    {
      for(int j = 0; j < Nx; ++j)
	{
	  temp[i * Nx + j] = T0;
	  
	  if(i == 0)
	    temp[i * Nx + j] = T1;
	  if(j == 0 && i != 0)
	    temp[i * Nx + j] = T2;
	  if(j == Nx-1)
	    temp[i * Nx + j] = T4;
	  if(i == Ny-1)
	    temp[i * Nx + j] = T3;
	}
    }

  // Fills positions
  pos.x.push_back(0.0);
  pos.y.push_back(0.0);
  for(int i = 1; i < Ny; ++i)
    {
      if(i < Nx)
	pos.x.push_back(pos.x[i-1] + dx);
      pos.y.push_back(pos.y[i-1] + dy);
    }

  // Fills the matrix
  for(int i = 0; i < Ax.n; ++i)
    {
      double alphax = diff/(dx*dx);
      double alphay = diff/(dy*dy);
      Ax.diag[i] = 1/dt_2 + 2*alphax;
      Ay.diag[i] = 1/dt_2 + 2*alphay;
      if(i < n-1)
	{
	  Ax.upper[i] = -alphax;
	  Ay.upper[i] = -alphay;
	}
      if(i > 0)
	{
	  Ax.lower[i] = -alphax;
	  Ay.lower[i] = -alphay;
	}
    }

  // Forces B.C.'s
  for(int i = 0; i < Ny; ++i)
    {
      for(int j = 0; j < Nx; ++j)
	{
	  if(i == 0)
	    {
	      Ax.diag[i * Nx + j] = 1;
	      Ax.upper[i * Nx + j] = 0;
	      Ax.lower[i * Nx + j] = 0;
	      Ay.diag[i * Nx + j] = 1;
	      Ay.upper[i * Nx + j] = 0;
	      Ay.lower[i * Nx + j] = 0;
	    }
	  if(j == 0 && i != 0)
	    {
	      Ax.diag[i * Nx + j] = 1;
	      Ax.upper[i * Nx + j] = 0;
	      Ax.lower[i * Nx + j] = 0;
	      Ay.diag[i * Nx + j] = 1;
	      Ay.upper[i * Nx + j] = 0;
	      Ay.lower[i * Nx + j] = 0;
	    }

	  if(j == Nx-1)
	    {
	      Ax.diag[i * Nx + j] = 1;
	      Ax.upper[i * Nx + j] = 0;
	      Ax.lower[i * Nx + j] = 0;
	      Ay.diag[i * Nx + j] = 1;
	      Ay.upper[i * Nx + j] = 0;
	      Ay.lower[i * Nx + j] = 0;
	    }

	  if(i == Ny-1)
	    {
	      Ax.diag[i * Nx + j] = 1;
	      Ax.upper[i * Nx + j] = 0;
	      Ax.lower[i * Nx + j] = 0;
	      Ay.diag[i * Nx + j] = 1;
	      Ay.upper[i * Nx + j] = 0;
	      Ay.lower[i * Nx + j] = 0;
	    }

	}
    }
  Ax.diag[0] = 1;
  Ay.diag[0] = 1;
  Ax.upper[0] = 0;
  Ay.upper[0] = 0;
  Ax.lower[0] = 0;
  Ay.lower[0] = 0;

  WriteFile(pos, temp, to_string(time));
 
  for(int i = 0; i <= num_itr; ++i)
    {
      time += dt;
      double alphax = diff/(dy*dy);
      double alphay = diff/(dx*dx);      
      // Build b for n+1/2
      for(int i = 1; i < Ny-1; ++i)
	{
	  for(int j = 1; j < Nx-1; ++j)
	    b[i * Nx + j] = (1/dt_2 - 2*alphax)*temp[i * Nx + j] + alphax * temp[i * Nx + (j+1)] + alphax * temp[i * Nx + (j-1)];
	}
      
      // Impose B.C.
      for(int i = 0; i < Ny; ++i)
	{
	  for(int j = 0; j < Nx; ++j)
	    {
	      if(i == 0)
		b[i * Nx + j] = T1;
	      else if(j == 0)
		b[i * Nx + j] = T2;
	      else if(j == Nx-1)
		b[i * Nx + j] = T4;
	      else if(i == Ny-1)
		b[i * Nx + j] = T3;
	    }
	}
      
      TridiagSolve(Ax, b, temp);

       // Build b for n+1
      for(int i = 1; i < Ny-1; ++i)
	{
	  for(int j = 1; j < Nx-1; ++j)
	    b[i * Nx + j] = (1/dt_2 - 2*alphay)*temp[i * Nx + j] + alphay * temp[(i+1) * Nx + j] + alphay * temp[(i-1) * Nx + j];
	}
      
      // Impose B.C.
      for(int i = 0; i < Ny; ++i)
	{
	  for(int j = 0; j < Nx; ++j)
	    {
	      if(i == 0)
		b[i * Nx + j] = T1;
	      else if(j == 0)
		b[i * Nx + j] = T2;
	      else if(j == Nx-1)
		b[i * Nx + j] = T4;
	      else if(i == Ny-1)
		b[i * Nx + j] = T3;
	    }
	}

      TridiagSolve(Ay, b, temp);
      
      if(WriteTime(time, writeTimes))
	WriteFile(pos, temp, to_string(time));
      
    }

  cout << "Searching for steady state" << endl;

  while(SteadyState(pos, old_temp, temp, steadyState) == false)
    {
      old_temp = temp;
      time += dt;
      double alphax = diff/(dy*dy);
      double alphay = diff/(dx*dx);      
      // Build b for n+1/2
      for(int i = 1; i < Ny-1; ++i)
	{
	  for(int j = 1; j < Nx-1; ++j)
	    b[i * Nx + j] = (1/dt_2 - 2*alphax)*temp[i * Nx + j] + alphax * temp[i * Nx + (j+1)] + alphax * temp[i * Nx + (j-1)];
	}
      
      // Impose B.C.
      for(int i = 0; i < Ny; ++i)
	{
	  for(int j = 0; j < Nx; ++j)
	    {
	      if(i == 0)
		b[i * Nx + j] = T1;
	      else if(j == 0)
		b[i * Nx + j] = T2;
	      else if(j == Nx-1)
		b[i * Nx + j] = T4;
	      else if(i == Ny-1)
		b[i * Nx + j] = T3;
	    }
	}
      
      TridiagSolve(Ax, b, temp);

       // Build b for n+1
      for(int i = 1; i < Ny-1; ++i)
	{
	  for(int j = 1; j < Nx-1; ++j)
	    b[i * Nx + j] = (1/dt_2 - 2*alphay)*temp[i * Nx + j] + alphay * temp[(i+1) * Nx + j] + alphay * temp[(i-1) * Nx + j];
	}
      
      // Impose B.C.
      for(int i = 0; i < Ny; ++i)
	{
	  for(int j = 0; j < Nx; ++j)
	    {
	      if(i == 0)
		b[i * Nx + j] = T1;
	      else if(j == 0)
		b[i * Nx + j] = T2;
	      else if(j == Nx-1)
		b[i * Nx + j] = T4;
	      else if(i == Ny-1)
		b[i * Nx + j] = T3;
	    }
	}

      TridiagSolve(Ay, b, temp);
    }

  WriteFile(pos, temp, to_string(time));
 
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

template <typename T>
void WriteFile(Dimensions<T> pos, vector<double> u, string time)
{
  stringstream ss;
  ss << time[0] << time[1] << time[2] << time[3] << time[4];
  string true_time = ss.str();
  cout << true_time << " seconds" << endl;
  string outputfile = true_time + "_sec.dat";
  string sys_cmd = "./plot_3.py " + true_time;
  ofstream outfile;
  outfile.open(outputfile.c_str());
  outfile << setprecision(6) << fixed;

  for(int i = 0; i < pos.y.size(); ++i)
    {
      for(int j = 0; j < pos.x.size(); ++j)
	outfile << pos.x[j] << ',' << pos.y[i] << ',' << u[i * pos.x.size() + j] << endl;
    }
  
  system(sys_cmd.c_str());

  return;
}

bool WriteTime(double time, vector<double> writeTimes, double tol)
{
  for(int i = 0; i < writeTimes.size(); ++i)
    {
      if(abs(time - writeTimes[i]) < tol)
	return true;
    }
  return false;
}

template <typename T>
bool SteadyState(Dimensions<T> pos, vector<double> T1, vector<double> T2, double steadyState)
{
  double sum = 0.0;
  for(int i = 0; i < T1.size(); ++i)
    {
      sum += abs(T1[i] - T2[i]);
    }
  sum = sum/(pos.x.size() * pos.y.size());
  if(sum < steadyState)
    return true;
  else
    return false;
}
