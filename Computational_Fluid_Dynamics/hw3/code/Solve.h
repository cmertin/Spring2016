#ifndef SOLVE_H
#define SOLVE_H

#include <iostream>
#include <vector>
#include <cmath>

template <typename T>
double Norm(std::vector<T> &oldVec, std::vector<T> &newVec)
{
  double err = 0.0;
  for(int i = 0; i < oldVec.size(); ++i)
    err += pow(abs(oldVec[i] - newVec[i]),2);

  return sqrt(err);
}

template <typename T>
void Jacobi(Flow<T> &flow, Dimensions<T> &dim, const int MAXITR=10000, const double &tol=10e-6)
{
  double error = 100;
  double dx2 = dim.dx * dim.dx;
  int Nx = dim.Nx;
  int Ny = dim.Ny;
  std::vector<T> psiOld = flow.psi;
  for(int itr = 0; itr < MAXITR; ++itr)
    {
      psiOld = flow.psi;
      for(int j = 1; j < Ny-1; ++j)
	{
	  for(int i = 1; i < Nx-1; ++i)
	    {
	      flow.psi[i + j*Nx] = 0.25 * (psiOld[i+1 + j*Nx] + psiOld[i-1 + j*Nx] + psiOld[i + (j+1)*Nx] + psiOld[i + (j-1)*Nx] + dx2 * flow.omega[i + j*Nx]);
	    }
	}
      if(itr % 10 == 0)
	error = Norm(flow.psi, psiOld);

      if(error < tol)
	return;
    }
  std::cerr << "JACOBI EXCEEDED MAXIMUM NUMBER OF ITERATIONS (" << MAXITR << ")" << std::endl;
  return;
}

template <typename T>
void GaussSeidel(Flow<T> &flow, Dimensions<T> &dim, const int MAXITR=10000, const double &tol = 10e-6)
{
  double error = 100;
  double dx2 = dim.dx * dim.dx;
  int Nx = dim.Nx;
  int Ny = dim.Ny;
  std::vector<T> psiOld = flow.psi;
  for(int itr = 0; itr < MAXITR; ++itr)
    {
      psiOld = flow.psi;
      for(int j = 1; j < Ny-1; ++j)
	{
	  for(int i = 1; i < Nx-1; ++i)
	    {
	      flow.psi[i + j*Nx] = 0.25 * (psiOld[i+1 + j*Nx] + flow.psi[i-1 + j*Nx] + psiOld[i + (j+1)*Nx] + flow.psi[i + (j-1)*Nx] + dx2 * flow.omega[i + j*Nx]);
	    }
	}
      if(itr % 10 == 0)
	error = Norm(flow.psi, psiOld);

      if(error < tol)
	return;
    }
  std::cerr << "GAUSS SEIDEL EXCEEDED MAXIMUM NUMBER OF ITERATIONS (" << MAXITR << ")" << std::endl;
  return;
}

template <typename T>
void GaussSeidel_SOR(Flow<T> &flow, Dimensions<T> &dim, double alpha, const int MAXITR=10000, const double &tol = 10e-6)
{
  double error = 100;
  double dx2 = dim.dx * dim.dx;
  int Nx = dim.Nx;
  int Ny = dim.Ny;
  double coeff = alpha * 0.25;
  std::vector<T> psiOld = flow.psi;
  for(int itr = 0; itr < MAXITR; ++itr)
    {
      psiOld = flow.psi;
      for(int j = 1; j < Ny-1; ++j)
	{
	  for(int i = 1; i < Nx-1; ++i)
	    {
	      flow.psi[i + j*Nx] = (1-alpha)*psiOld[i + j*Nx] + coeff * (psiOld[i+1 + j*Nx] + flow.psi[i-1 + j*Nx] + psiOld[i + (j+1)*Nx] + flow.psi[i + (j-1)*Nx] + dx2 * flow.omega[i + j*Nx]);
	    }
	}
      if(itr % 10 == 0)
	error = Norm(flow.psi, psiOld);

      if(error < tol)
	return;
    }
  std::cerr << "GAUSS SEIDEL [SOR] EXCEEDED MAXIMUM NUMBER OF ITERATIONS (" << MAXITR << ")" << std::endl;
  return;
}

template <typename T>
void MultiGrid(Flow<T> &flow, Dimensions<T> &dim)
{
  return;
}

#endif
