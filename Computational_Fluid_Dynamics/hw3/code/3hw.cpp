#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>
#include <string>
#include <algorithm>
#include "Dimensions.h"
#include "Output.h"
#include "Matrix.h"
//#include "Multigrid.h"

using namespace std;

template <typename T>
void UpdateVorticity(vector<T> &omega, vector<T> &u, vector<T> &v, const int &Nx, const double &U, const double &h);
template <typename T>
void VorticityWall(vector<T> &omega, vector<T> &psi, double &U, double &dx, int &Nx);
template <typename T>
void VorticityTime(vector<T> &omega, vector<T> &omegaOld, vector<T> &psi, double dx, double nu, double dt, int Nx);
template <typename T>
void UpdateVelocities(vector<T> &psi, vector<T> &u, vector<T> &v, double &dx, double &U, int &Nx);
template <typename T>
void VorticityPsi(vector<T> &omega, vector<T> &psi, int &Nx, double &dx);
template <typename T>
bool SteadyState(vector<T> &omega, vector<T> &omegaOld, int &Nx, int &Ny, double &tol=10e-10);

int main(int argc, char *argv[])
{
  Dimensions<double> pos;
  double Re = 100;
  double nu = 0.1;
  double minX = 0.0;
  double maxX = 1.0;
  double minY = 0.0;
  double maxY = 1.0;
  double U = (Re*nu)/(maxX - minX);
  const int MAXITR = 1000;
  int Nx = 10;//50;
  int Ny = 10;//50;
  double dx = (maxX - minX)/(Nx - 1);
  double dy = (maxY - minY)/(Ny - 1);
  double tol = 10e-10;
  double dt = 0.01;//dx/(40*U);
  vector<double> u(Nx * Ny, 0);
  vector<double> v(Nx * Ny, 0);
  vector<double> omega(Nx * Ny, 0);
  vector<double> omegaOld(Nx * Ny, 0);
  vector<double> psi(Nx * Ny, 0);
  vector<double> psiOld(Nx * Ny, 0);
  vector<double> psiOldCheck(Nx * Ny, 0);
  vector<double> residual(Nx * Ny, 0);
  vector<double> residualOld(Nx * Ny, 0);
  Matrix<double> A;
  A.diag.resize(Nx * Nx, 0);
  A.upper.resize(Nx * Nx, 0);
  A.uupper.resize(Nx * Nx, 0);
  A.lower.resize(Nx * Nx, 0);
  A.llower.resize(Nx * Nx, 0);
  A.Nx = Nx;
  A.Ny = Ny;
  double time = 0.0;
  int counter = 0;
  double error = 100;
  bool done = false;

  // Builds the x,y "grid"
  pos.x.resize(Nx);
  pos.y.resize(Ny);
  for(int i = 1; i < pos.x.size(); ++i)
    {
      pos.x[i] = pos.x[i-1] + dx;
      pos.y[i] = pos.y[i-1] + dy;
    }

  // Sets the velocity on the grid for the movement
  for(int i = 0; i < u.size(); ++i)
    {
      if(i > u.size() - Nx - 1)
	u[i] = U;
    }

  // Build Initial A Matrix
  for(int j = 0; j < Ny; ++j)
    {
      for(int i = 0; i < Nx; ++i)
	{
	  int index = i + j * Nx;
	  A.diag[index] = 4;
	  if(index >= Nx)
	    A.llower[index] = -1; 
	  if(index < Nx*Nx - Nx - 1)
	    A.uupper[index] = -1;
	  if(index < Nx*Nx - 1)
	    A.upper[index] = -1; 
	  if(index > 0)
	    A.lower[index] = -1;
	  if(i == Nx - 1)
	    A.upper[index] = 0;
	  if(i == 0 && j > 0)
	    A.lower[index] = 0;
	}
    }

    // Builds initial omega
    UpdateVorticity(omega, u, v, Nx, U, dx);

    // Get psi for T = 0
    for(int itr = 0; itr < MAXITR; ++itr)
       {
	 psiOld = psi;
	 for(int i = 0; i < Nx; ++i)
	   {
	     double dx2 = dx*dx;
	     for(int j = 0; j < Nx; ++j)
	       {
		 if(i == 0 || j == 0 || i == Nx-1 || j == Nx-1)
		   psi[j + i*Nx] = 0;
		 else
		   psi[j + i*Nx] = 0.25 * (psiOld[j+1 + i*Nx] + psiOld[j-1 + i*Nx] + psiOld[j + (i+1)*Nx] + psiOld[j + (i-1)*Nx] + dx2*omega[j + i*Nx]);
	       }
	   }
       }

    VorticityWall(omega, psi, U, dx, Nx);

    ++counter;

    // Iterates through time
   do{
     time += dt;
     cout << "Time: " << time << "\t(" << counter+1 << ')' <<  endl;

     psiOldCheck = psi;
     omegaOld = omega;
     
     // Update inner vorticity in time
     VorticityTime(omega, omegaOld, psi, U, dx, nu, Nx);
     
     // Solve for Psi with new values of omega
     // Jacobi(A, psi, omega, MAXITR);
     for(int itr = 0; itr < MAXITR; ++itr)
       {
	 double dx2 = dx*dx;
	 psiOld = psi;
	 for(int i = 1; i < Nx-1; ++i)
	   {
	     for(int j = 1; j < Nx-1; ++j)
	       {
		 psi[j + i*Nx] = 0.25 * (psiOld[j+1 + i*Nx] + psiOld[j-1 + i*Nx] + psiOld[j + (i+1)*Nx] + psiOld[j + (i-1)*Nx] + dx2*omega[j + i*Nx]);
	       }
	   }
	 if(Norm(psi, psiOld) < tol && itr % 10 == 0 && itr > 0)
	   {
	     cout << "ITR: " << itr << endl;
	     break;
	   }
       }

     // Find the velocity components from Psi
     UpdateVelocities(psi, u, v, dx, U, Nx);

     // Update the vorticity at the wall and interior
     VorticityPsi(omega, psi, Nx, dx);
     VorticityWall(omega, psi, U, dx, Nx);
     
     // Calculate the error
     error = Norm(omega, omegaOld);
     error = SteadyState(omega, omegaOld, Nx, Ny, tol);

     cout << "Error: " << error << endl;
     
     ++counter;
   }while(error > tol);

  string fileout = "finished.dat";
  WriteFile(pos, omega, psi, u, v, fileout);

  cout << "Plotting... ";
  cout.flush();
  system("./plot.py");
  cout << "DONE" << endl;
  
  return 0;
}

template <typename T>
void UpdateVorticity(vector<T> &omega, vector<T> &u, vector<T> &v, const int &Nx, const double &U, const double &h)
{
  double h_ = h*2;
  for(int i = 0; i < Nx; ++i)
    {
      for(int j = 0; j < Nx; ++j)
	{
	  double du_dy;
	  double dv_dx;
	  if(i == 0) // Bottom boundary
	    {
	      du_dy = 0;
	      dv_dx = 0;
	    }
	  else if(i == Nx - 1) // Top boundary
	    {
	      du_dy = U;
	      dv_dx = 0;
	    }
	  else if(j == 0) // Left boundary
	    {
	      du_dy = 0;
	      dv_dx = 0;
	    }
	  else if(j == Nx-1) // Right boundary
	    {
	      du_dy = 0;
	      dv_dx = 0;
	    }
	  else // Inner
	    {
	      du_dy = (u[j + (i+1)*Nx] - u[j + (i-1)*Nx])/h_;
	      dv_dx = (v[j+1 + i*Nx] - v[j-1 + i*Nx])/h_;
	    }
	  omega[j + i * Nx] = dv_dx - du_dy;
	}
    }
  return;
}

template <typename T>
void VorticityWall(vector<T> &omega, vector<T> &psi, double &U, double &dx, int &Nx)
{
  double dx2 = dx*dx;
  double coeff = -(2.0/dx2);
  for(int i = 0; i < Nx; ++i)
    {
      for(int j = 0; j < Nx; ++j)
	{
	  if(i == Nx - 1) // Top boundary
	    omega[j + i*Nx] = coeff*(psi[j + (i-1)*Nx] + U*dx);
	  else if(i == 0) // Bottom boundary
	    omega[j + i*Nx] = coeff*(psi[j + (i+1)*Nx]);
	  else if(j == 0) // Left boundary
	    omega[j + i*Nx] = coeff*(psi[j+1 + i*Nx]);
	  else if(j == Nx - 1) // Right boundary
	    omega[j + i*Nx] = coeff*(psi[j-1 + i*Nx]);
	}
    }
  return;
}

template <typename T>
void VorticityTime(vector<T> &omega, vector<T> &omegaOld, vector<T> &psi, double dx, double nu, double dt, int Nx)
{
  double dx_ = 2.0*dx;
  double dx2 = dx*dx;
  for(int i = 1; i < Nx-1; ++i)
    {
      for(int j = 1; j < Nx-1; ++j)
	{
	  double dpsi_dx = (psi[j+1 + i*Nx] - psi[j-1 + i*Nx])/(dx_);
	  double dpsi_dy = (psi[j + (i+1)*Nx] - psi[j + (i-1)*Nx])/(dx_);
	  double domega_dx = (omegaOld[j+1 + i*Nx] - omegaOld[j-1 + i*Nx])/(dx_);
	  double domega_dy = (omegaOld[j + (i+1)*Nx] - omegaOld[j + (i-1)*Nx])/(dx_);
	  double omega_old = omegaOld[j + i*Nx];
	  double omega_dxp = omegaOld[j+1 + i*Nx];
	  double omega_dxm = omegaOld[j-1 + i*Nx];
	  double omega_dyp = omegaOld[j + (i+1)*Nx];
	  double omega_dym = omegaOld[j + (i-1)*Nx];
	  double omega_d2 = (omega_dxp + omega_dxm + omega_dyp + omega_dym - 4*omega_old)/dx2;
	  omega[j + i*Nx] = omega_old + dt * (dpsi_dy * domega_dx - dpsi_dx * domega_dy + nu * omega_d2);
	}
    }
  return;
}

template <typename T>
void UpdateVelocities(vector<T> &psi, vector<T> &u, vector<T> &v, double &dx, double &U, int &Nx)
{
  double dx2 = dx*dx;
  for(int i = 0; i < Nx; ++i)
    {
      for(int j = 0; j < Nx; ++j)
	{
	  int index = j + i*Nx;
	  double dpsi_dx = (psi[j+1 + i*Nx] - 2*psi[index] + psi[j-1 + i*Nx])/dx2;
	  double dpsi_dy = (psi[j + (i+1)*Nx] - 2*psi[index] + psi[j + (i-1)*Nx])/dx2;
	  if(i == 0) // Bottom Boundary
	    {
	      u[index] = 0;
	      v[index] = 0;
	    }
	  else if(i == Nx-1) // Top Boundary
	    {
	      u[index] = U;
	      v[index] = 0;
	    }	      
	  else if(j == 0) // Left Boundary
	    {
	      u[index] = 0;
	      v[index] = 0;
	    }
	  else if(j == Nx-1) // Right Boundary
	    {
	      u[index] = 0;
	      v[index] = 0;
	    }
	  else // Inside
	    {
	      u[index] = dpsi_dy;
	      v[index] = -dpsi_dx;
	    }
	}
    }
  return;
}

template <typename T>
void VorticityPsi(vector<T> &omega, vector<T> &psi, int &Nx, double &dx)
{
  double dx2 = dx*dx;
  for(int i = 1; i < Nx-1; ++i)
    {
      for(int j = 1; j < Nx-1; ++j)
	{
	  double d2psi_dx2 = (psi[j+1 + i*Nx] - 2 * psi[j + i*Nx] + psi[j-1 + i*Nx])/dx2;
	  double d2psi_dy2 = (psi[j + (i+1)*Nx] - 2*psi[j + i*Nx] + psi[j + (i-1)*Nx])/dx2;
	  omega[j + i*Nx] = - (d2psi_dx2 + d2psi_dy2);
	}
    }
  return;
}

template <typename T>
bool SteadyState(vector<T> &omega, vector<T> &omegaOld, int &Nx, int &Ny, double &tol)
{
  T sum = 0.0;
  for(int i = 0; i < omega.size(); ++i)
    {
      sum += abs(omega[i] - omegaOld[i]);
    }
  //sum = 1/(Nx * Ny) * sum;
  if(sum < tol)
    return true;
  else
    return false;
}
