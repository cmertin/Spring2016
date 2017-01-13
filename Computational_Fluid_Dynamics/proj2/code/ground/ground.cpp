#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include <omp.h>
#include <ctime>
#include <chrono>
#include "Dimensions.h"

using namespace std;

template <typename T>
void BoundaryConditions(Flow<T> &flow, Dimensions<T> &dim, Train<T> &train, int threads=1);
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
  if(argc != 7)
    {
      cerr << "\nNote: The following arguments apply to the train" << endl;
      cerr << "Usage: ./proj2 U_train startX startY width height separation\n" << endl;
      return -1;
    }
  
  // Initialize variables
  chrono::time_point<chrono::system_clock> start,end;
  start = chrono::system_clock::now();
  Flow<double> flow;
  Dimensions<double> dim;
  double time = 0.0;
  double alpha = 1.25;
  double tol = 10e-6;
  double error = 100;
  int threads = omp_get_max_threads();
  const int MAX_ITR = 100000;
  int counter = -1;
  string filename = "train_";
  string filename_arg = "";
  string sys_cmd = "./plot.py ";
  // Train specifications
  Train<double> train;
  train.velocity = atoi(argv[1]);
  train.startX = atoi(argv[2]);
  train.startY = atoi(argv[3]);
  train.width = atoi(argv[4]);
  train.height = atoi(argv[5]);
  train.separation = atoi(argv[6]);  
  // Initialize Dimensions
  dim.Nx = 100;
  dim.Ny = 100;
  dim.maxX = 2.0;
  dim.minX = 0.0;
  dim.maxY = 2.0;
  dim.minY = 0.0;
  dim.dx = (dim.maxX - dim.minX)/(dim.Nx - 1);
  dim.dy = (dim.maxY - dim.minY)/(dim.Ny - 1);
  dim.x.resize(dim.Nx, 0);
  dim.y.resize(dim.Ny, 0);
  dim.Re = 100;
  // Initialize Flows
  flow.u.resize(dim.Nx+1, vector<double>(dim.Ny+2, 0));
  flow.ut.resize(dim.Nx+1, vector<double>(dim.Ny+2, 0));
  flow.v.resize(dim.Nx+2, vector<double>(dim.Ny+1, 0));
  flow.vt.resize(dim.Nx+2, vector<double>(dim.Ny+1, 0));
  flow.P.resize(dim.Nx+2, vector<double>(dim.Ny+2, 0));
  vector<vector<double> > Pold = flow.P;
  flow.omega.resize(dim.Nx, vector<double>(dim.Ny, 0));
  flow.U_top = 0;
  flow.U_bottom = 0;
  flow.U_left = train.velocity;
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
      ++counter;
      cout << setprecision(4) << fixed << "Time: " << time
	   << "\t(" << counter << ')' << endl;

      if(counter % 100 == 0)
	Pold = flow.P;

      // Impose Boundary Conditions for Velocity using Reflexive
      BoundaryConditions(flow, dim, train, threads);

      TemporaryVelocities(flow, dim, threads);

      // Solve for Pressure
      GaussSeidel_SOR(flow, dim, alpha, MAX_ITR);
 
      CorrectVelocities(flow, dim, threads);

      if(counter % 100 == 0)
	error = Norm(flow.P, Pold);

      time += dim.dt;
    }while(error > tol);

  BoundaryConditions(flow, dim, train, threads);

  GetVorticity(flow, dim, threads);

  // Output the data
  end = chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;
  time_t start_time = chrono::system_clock::to_time_t(start);
  time_t end_time = chrono::system_clock::to_time_t(end);
  filename.append(to_string((int)train.velocity));
  filename.append("_");
  filename.append(to_string(train.width));
  filename.append("_");
  filename.append(to_string(train.height));
  filename.append("_");
  filename.append(to_string(train.separation));
  filename_arg = filename;
  filename_arg.append("_data");
  filename.append("_data.dat");
  cout << "Output File: " << filename << endl;
  ofstream outfile;
  outfile.open(filename.c_str());
  outfile << setprecision(6) << fixed;
  outfile << "# ================" << endl;
  outfile << "# Train Parameters" << endl;
  outfile << "# ================" << endl;
  outfile << "# Train Velocity: " << (int)train.velocity << endl;
  outfile << "# Start X: " << train.startX << endl;
  outfile << "# Start Y: " << train.startY << endl;
  outfile << "# Width: " << train.width << endl;
  outfile << "# Height: " << train.height << endl;
  outfile << "# Separation: " << train.separation << endl;
  outfile << "# =================" << endl;
  outfile << "# Domain Parameters" << endl;
  outfile << "# =================" << endl;
  outfile << "# min X, max X: " << dim.minX << ',' << dim.maxX << endl;
  outfile << "# min Y, max Y: " << dim.minY << ',' << dim.maxY << endl;
  outfile << "# Nx: " << dim.Nx << endl;
  outfile << "# Ny: " << dim.Ny << endl;
  outfile << "# dx: " << dim.dx << endl;
  outfile << "# dy: " << dim.dy << endl;
  outfile << "# dt: " << dim.dt << endl;
  outfile << "# Re: " << dim.Re << endl;
  outfile << "# nu: " << dim.nu << endl;
  outfile << "# U: " << dim.U << endl;
  outfile << "# Final time: " << time << endl;
  outfile << "# Iterations: " << counter << endl;
  outfile << "# Start Time: " << ctime(&start_time);
  outfile << "# End Time: " << ctime(&end_time);
  outfile << "# Duration: " << elapsed_seconds.count() << " seconds" << endl;
  outfile << "# ====" << endl;
  outfile << "# Data" << endl;
  outfile << "# ====" << endl;
  outfile << "# x\t y\t u\t v\t P\t vorticity" << endl;
  for(int j = 0; j < dim.x.size(); ++j)
    {
      for(int i = 0; i < dim.y.size(); ++i)
	{
	  outfile << dim.x[i] << ',' << dim.y[j] << ',' << flow.u[i+1][j+1] << ',' << flow.v[i+1][j] << ',' << flow.P[i+1][j+1]  << ',' << flow.omega[i][j] << endl;
	}
    }
  outfile.close();

  sys_cmd.append(filename_arg);
  sys_cmd.append(" ");
  sys_cmd.append(to_string((int)train.velocity));
  sys_cmd.append(" ");
  sys_cmd.append(to_string(train.startX));
  sys_cmd.append(" ");
  sys_cmd.append(to_string(train.startY));
  sys_cmd.append(" ");
  sys_cmd.append(to_string(train.width));
  sys_cmd.append(" ");
  sys_cmd.append(to_string(train.height));
  sys_cmd.append(" ");
  sys_cmd.append(to_string(dim.dx));
  sys_cmd.append(" ");
  sys_cmd.append(to_string(train.separation));

  cout << "Plotting... ";
  cout.flush();
  system(sys_cmd.c_str());
  cout << "DONE" << endl;
  
  return 0;
}

template <typename T>
void BoundaryConditions(Flow<T> &flow, Dimensions<T> &dim, Train<T> &train, int threads)
{
  /* Inlet
     =====
     u = U
     v = 0
     P = 0
  */
  #pragma omp parallel num_threads(threads)
  for(int j = 1; j < dim.Ny+1; ++j)
    {
      flow.u[0][j] = flow.U_left;
      flow.v[0][j] = 0;
      flow.P[0][j] = -flow.P[1][j];
    }

  /* Outlet
     ======
     du/dx = 0 (backward difference)
     dv/dx = 0 (backward difference)
     dP/dx = 0 (backward difference)
  */
  #pragma omp parallel num_threads(threads)
  for(int j = 1; j < dim.Ny+1; ++j)
    {
      flow.u[dim.Nx][j] = -flow.u[dim.Nx-1][j];
      flow.v[dim.Nx+1][j] = -flow.v[dim.Nx][j];
      flow.P[dim.Nx+1][j] = -flow.P[dim.Nx][j];
    }

  /* Periodic B.C.
     =============
     P(top,ghost) = P(bottom)
     P(bottom,ghost) = P(top)
     v(bottom) = v(top)
     u(top,ghost) = u(bottom)
     u(bottom,ghost) = u(top)
  */
  #pragma omp parallel num_threads(threads)
  for(int i = 0; i < dim.Nx+2; ++i)
    {
      flow.P[i][dim.Ny+1] = flow.P[i][1];
      flow.P[i][0] = flow.P[i][dim.Ny];
      flow.v[i][0] = flow.v[i][dim.Ny];
      if(i < dim.Nx+1)
	{
	  flow.u[i][dim.Ny+1] = flow.u[i][1];
	  flow.u[i][0] = flow.u[i][dim.Ny];
	}
    }

  // B.C. for train
  #pragma omp parallel num_threads(threads)
  for(int i = train.startX; i < train.startX + train.width; ++i)
    {
      // i for second car
      int i2 = i + train.width + train.separation;
      for(int j = train.startY; j < train.startY + train.height; ++j)
	{
	  // Front car
          flow.u[i][j] = 0;
	  flow.v[i][j] = 0;
          flow.P[i][j] = 0;
	  // Second car
	  //flow.u[i2][j] = 0;
	  //flow.v[i2][j] = 0;
	  //flow.P[i2][j] = 0;
	  // left boundary
	  if(i == train.startX)
	    {
	      flow.u[i][j] = -flow.u[i-1][j];
	      //flow.u[i2][j] = -flow.u[i2-1][j];
	    }
	  // right boundary
	  if(i == train.startX + train.width - 1)
	    {
	      flow.u[i][j] = flow.u[i+1][j];
	      //flow.u[i2][j] = -flow.u[i2+1][j];
	    }
	  // bottom boundary
	  if(j == train.startY)
	    {
	      flow.v[i][j] = -flow.v[i][j-1];
	      //flow.v[i2][j] = -flow.v[i2][j-1];
	    }
	  // top boundary
	  if(j == train.startY + train.height - 1)
	    {
	      flow.v[i][j] = flow.v[i][j+1];
	      //flow.v[i2][j] = flow.v[i2][j+1];
	    }
	}
    }

  // B.C. for ground
  for(int j = train.startY-5; j < train.startY+1; ++j)
    {
      for(int i = 0; i < dim.Nx; ++i)
	{
	  flow.u[i][j] = 0;
	  flow.v[i][j] = 0;
	  flow.P[i][j] = 0;
	  
	  if(j == train.startY-5)
	      flow.v[i][j] = -flow.v[i][j-1];

	  if(j == train.startY)
	    flow.v[i][j] = flow.v[i][j+1];
	}      
    }
      
    /**/
  /*
  // Build Boundary Conditions
  #pragma omp parallel num_threads(threads)
  for(int i = 0; i < dim.Nx+1; ++i)
    {
      flow.u[i][0] = 2 * flow.U_bottom - flow.u[i][1];
      flow.u[i][dim.Ny+1] = 2 * flow.U_top - flow.u[i][dim.Ny];
    }
  /*
  #pragma omp parallel num_threads(threads)
  for(int i = 0; i < dim.Ny+1; ++i)
    {
      flow.v[0][i] = 2 * flow.U_left - flow.v[1][i];
      flow.v[dim.Nx+1][i] = 2 * flow.U_right - flow.v[dim.Nx][i];
    }
  /**/
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
