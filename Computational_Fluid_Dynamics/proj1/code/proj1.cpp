#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <omp.h>
#include "Dimensions.h"

using namespace std;

template <typename T>
void BoundaryConditions(Flow<T> &flow, Dimensions<T> &dim, int threads=1);
template <typename T>
void TemporaryVelocities(Flow<T> &flow, Dimensions<T> &dim, int threads=1);
template <typename T>
T Norm(vector<vector<T> > &newVals, vector<vector<T> > &oldVals);
template <typename T>
void GaussSeidel_SOR(Flow<T> &flow, Dimensions<T> &dim, double alpha, const int MAX_ITR=10000, const double &tol = 10e-6);
template <typename T>
void CorrectVelocities(Flow<T> &flow, Dimensions<T> &dim, int threads=1);
template <typename T>
void GetVorticity(Flow<T> &flow, Dimensions<T> &dim, int threads=1);

int main(int argc, char *argv[])
{
  // Initialize variables
  Flow<double> flow;
  Dimensions<double> dim;
  double time = 0.0;
  double alpha = 1.2;
  double tol = 10e-6;
  double error = 100;
  int threads = omp_get_max_threads();
  const int MAX_ITR = 100000;
  int counter = 0;
  string filename = "Re=";
  string sys_cmd = "./plot.py ";
  // Initialize Dimensions
  dim.Nx = 50;
  dim.Ny = 50;
  dim.maxX = 1.0;
  dim.minX = 0.0;
  dim.maxY = 1.0;
  dim.minY = 0.0;
  dim.dx = (dim.maxX - dim.minX)/(dim.Nx - 1);
  dim.dy = (dim.maxY - dim.minY)/(dim.Ny - 1);
  dim.x.resize(dim.Nx, 0);
  dim.y.resize(dim.Ny, 0);
  if(argc == 1)
    dim.Re = 100;
  else
    dim.Re = atoi(argv[1]);

  // Initialize Flows
  flow.u.resize(dim.Nx+1, vector<double>(dim.Ny+2, 0));
  flow.ut.resize(dim.Nx+1, vector<double>(dim.Ny+2, 0));
  flow.v.resize(dim.Nx+2, vector<double>(dim.Ny+1, 0));
  flow.vt.resize(dim.Nx+2, vector<double>(dim.Ny+1, 0));
  flow.P.resize(dim.Nx+2, vector<double>(dim.Ny+2, 0));
  vector<vector<double> > Pold = flow.P;
  flow.omega.resize(dim.Nx, vector<double>(dim.Ny, 0));
  flow.U_top = 10;
  flow.U_bottom = 0;
  flow.U_left = 0;
  flow.U_right = 0;
  double vel = max(flow.U_top, flow.U_bottom) + max(flow.U_left, flow.U_right);
  dim.nu = vel/dim.Re;
  dim.dt = (dim.dx*dim.dx)/(40*dim.nu);


  // Build grid of x,y values
  for(int i = 1; i < dim.x.size(); ++i)
    dim.x[i] = dim.x[i-1] + dim.dx;

  for(int i = 1; i < dim.y.size(); ++i)
    dim.y[i] = dim.y[i-1] + dim.dy;


  do
    {
      cout << setprecision(4) << fixed << "Time: " << time
	   << "\t(" << counter << ')' << endl;

      if(counter % 100 == 0)
	Pold = flow.P;

      // Impose Boundary Conditions for Velocity using Reflexive
      BoundaryConditions(flow, dim, threads);

      TemporaryVelocities(flow, dim, threads);

      // Solve for Pressure
      GaussSeidel_SOR(flow, dim, alpha, MAX_ITR);
 
      CorrectVelocities(flow, dim, threads);

      if(counter % 100 == 0)
	error = Norm(flow.P, Pold);

      time += dim.dt;
      ++counter;
    }while(error > tol);

  GetVorticity(flow, dim, threads);

  // Output the data
  filename.append(to_string((int)dim.Re));
  filename.append("_data.dat");
  cout << "Output File: " << filename << endl;
  ofstream outfile;
  outfile.open(filename.c_str());
  outfile << setprecision(6) << fixed;
  for(int j = 0; j < dim.x.size(); ++j)
    {
      for(int i = 0; i < dim.y.size(); ++i)
	{
	  outfile << dim.x[i] << ',' << dim.y[j] << ',' << flow.u[i+1][j+1] << ',' << flow.v[i+1][j] << ',' << flow.P[i+1][j+1]  << ',' << flow.omega[i][j] << endl;
	}
    }
  outfile.close();

  sys_cmd.append(to_string((int)dim.Re));

  cout << "Plotting... ";
  cout.flush();
  system(sys_cmd.c_str());
  cout << "DONE" << endl;
  
  return 0;
}

template <typename T>
void BoundaryConditions(Flow<T> &flow, Dimensions<T> &dim, int threads)
{
  // Build Boundary Conditions
  #pragma omp parallel num_threads(threads)
  for(int i = 0; i < dim.Nx+1; ++i)
    {
      flow.u[i][0] = 2 * flow.U_bottom - flow.u[i][1];
      flow.u[i][dim.Ny+1] = 2 * flow.U_top - flow.u[i][dim.Ny];
    }
  
  #pragma omp parallel num_threads(threads)
  for(int i = 0; i < dim.Ny+1; ++i)
    {
      flow.v[0][i] = 2 * flow.U_left - flow.v[1][i];
      flow.v[dim.Nx+1][i] = 2 * flow.U_right - flow.v[dim.Nx][i];
    }
  return;
}

template <typename T>
void TemporaryVelocities(Flow<T> &flow, Dimensions<T> &dim, int threads)
{
  // Temporary u velocity
  #pragma omp parallel num_threads(threads)
  for(int i = 1; i < dim.Nx; ++i)
    {
      double coeff1 = -(0.25/dim.dx);
      double coeff2 = dim.nu/(dim.dx * dim.dx);
      for(int j = 1; j < dim.Ny+1; ++j)
	flow.ut[i][j] = flow.u[i][j] + dim.dt * (coeff1*(pow((flow.u[i+1][j]+flow.u[i][j]),2) - pow((flow.u[i][j] + flow.u[i-1][j]),2) + (flow.u[i][j+1] + flow.u[i][j])*(flow.v[i+1][j] + flow.v[i][j]) - (flow.u[i][j] + flow.u[i][j-1]) * (flow.v[i+1][j-1] + flow.v[i][j-1])) + coeff2 * (flow.u[i+1][j] + flow.u[i-1][j] + flow.u[i][j+1] + flow.u[i][j-1] - 4 * flow.u[i][j]));
    }
  
  // Temporary v velocity 
  #pragma omp paralle num_threads(threads)
  for(int i = 1; i < dim.Nx+1; ++i)
    {
      double coeff1 = -(0.25/dim.dx);
      double coeff2 = dim.nu/(dim.dx * dim.dx);
      for(int j = 1; j < dim.Ny; ++j)
	{
	  flow.vt[i][j] = flow.v[i][j] + dim.dt * (coeff1*((flow.u[i][j+1]+flow.u[i][j])*(flow.v[i+1][j] + flow.v[i][j]) - (flow.u[i-1][j+1] + flow.u[i-1][j]) * (flow.v[i][j] + flow.v[i-1][j]) + pow((flow.v[i][j+1] + flow.v[i][j]),2) - pow((flow.v[i][j] + flow.v[i][j-1]),2)) + coeff2 * (flow.v[i+1][j] + flow.v[i-1][j] + flow.v[i][j+1] + flow.v[i][j-1] - 4 * flow.v[i][j]));
	}
    }
  return;
}

template <typename T>
T Norm(std::vector<std::vector<T> > &newVals, std::vector<std::vector<T> > &oldVals)
{
  double err = 0.0;
  for(int i = 0; i < newVals.size(); ++i)
    {
      for(int j = 0; j < newVals[i].size(); ++j)
	err += pow(abs(newVals[i][j] - oldVals[i][j]),2);
    }

  return sqrt(err);
}

template <typename T>
void GaussSeidel_SOR(Flow<T> &flow, Dimensions<T> &dim, double alpha, const int MAX_ITR, const double &tol)
{
  std::vector<std::vector<double> > Pold = flow.P;
  double error = 100;
  double coeff = 0;
  double delta = (dim.dx/dim.dt);

  for(int itr = 0; itr < MAX_ITR; ++itr)
    {
      Pold = flow.P;
      for(int i = 1; i < dim.Nx+1; ++i)
	{
	  for(int j = 1; j < dim.Ny+1; ++j)
	    {
	      if((i == 1 && j == 1) || (i == dim.Nx && j == 1) ||
		 (i == 1 && j == dim.Ny) || (i == dim.Nx && j == dim.Ny))
		coeff = (1.0/2.0) * alpha;
	      else if(i == 1 || j == 1 || i == dim.Nx || j == dim.Ny)
		coeff = (1.0/3.0) * alpha;
	      else
		coeff = (1.0/4.0) * alpha;
	      flow.P[i][j] = (1 - alpha) * Pold[i][j] + coeff * (Pold[i+1][j] + flow.P[i-1][j] + Pold[i][j+1] + flow.P[i][j-1] - delta * (flow.ut[i][j] - flow.ut[i-1][j] + flow.vt[i][j] - flow.vt[i][j-1]));
	    }
	}
      if(itr % 10 == 0)
	error = Norm(flow.P, Pold);

      //std::cout << "Error: " << error << std::endl;
      if(error < tol)
	return;
    }
  std::cerr << "GAUSS SEIDEL [SOR] EXCEEDED MAXIMUM NUMBER OF ITERATIONS (" << MAX_ITR << ")" << std::endl;
  return;
}

template <typename T>
void CorrectVelocities(Flow<T> &flow, Dimensions<T> &dim, int threads)
{
  // Correct the velocities
  #pragma omp parallel num_threads(threads)
  for(int i = 1; i < dim.Nx; ++i)
    {
      double delta = dim.dt/dim.dx;
      for(int j = 1; j < dim.Ny+1; ++j)
	flow.u[i][j] = flow.ut[i][j] - delta * (flow.P[i+1][j] - flow.P[i][j]);
    }
  
  #pragma omp parallel num_threads(threads)
  for(int i = 1; i < dim.Nx+1; ++i)
    {
      double delta = dim.dt/dim.dx;
      for(int j = 1; j < dim.Ny; ++j)
	flow.v[i][j] = flow.vt[i][j] - delta * (flow.P[i][j+1] - flow.P[i][j]);
    }
  return;
}

template <typename T>
void GetVorticity(Flow<T> &flow, Dimensions<T> &dim, int threads)
{
  // Calculate vorticity
  #pragma omp parallel num_threads(threads)
  for(int i = 1; i < dim.Nx-1; ++i)
    {
      double dx_ = dim.dx * 2.0;
      double dy_ = dim.dy * 2.0;
      for(int j = 1; j < dim.Ny-1; ++j)
	{
	  double du_dy = (flow.u[i][j+1] - flow.u[i][j])/dy_;
	  double dv_dx = (flow.v[i+1][j] - flow.v[i][j])/dx_;
	  flow.omega[i][j] = dv_dx - du_dy;
	}
    }
  return;
}
