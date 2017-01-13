#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

template <typename T>
struct Matrix
{
  std::vector<T> diag;
  std::vector<T> upper;
  std::vector<T> lower;
  std::vector<T> uupper;
  std::vector<T> llower;
  int Nx;
  int Ny;
};

template <typename T>
double Norm(std::vector<T> &oldVec, std::vector<T> &newVec)
{
  double err = 0.0;
  for(int i = 0; i < oldVec.size(); ++i)
    err += pow(abs(oldVec[i] - newVec[i]),2);

  return sqrt(err);
}

template <typename T>
void MatVec(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b)
{
  assert(A.diag.size() == x.size());
  for(int i = 0; i < x.size(); ++i)
    {
      b[i] = A.diag[i] * x[i];
      if(i > 0)
	b[i] += A.lower[i] * x[i - 1];
      if(i >= A.Nx)
	b[i] += A.llower[i] * x[i - A.Nx];
      if(i < A.Nx * A.Nx)
	b[i] += A.upper[i] * x[i + 1];
      if(i < A.Nx * A.Nx - A.Nx)
	b[i] += A.upper[i] * x[i + A.Nx];
    }
  return;
}

template <typename T>
void VecVec(std::vector<T> &a, std::vector<T> &x, T &b)
{
  assert(a.size() == x.size());
  b = 0.0;
  for(int i = 0; i < x.size(); ++i)
    b += a[i] * x[i];
  
  return;
}

template <typename T>
void VecSubVec(std::vector<T> &a, std::vector<T> &x, std::vector<T> &b)
{
  assert(a.size() == x.size());
  for(int i = 0; i < x.size(); ++i)
    b[i] = a[i] - x[i];
    
  return;
}

template <typename T>
void VecPlusVec(std::vector<T> &a, std::vector<T> &x, std::vector<T> &b)
{
  assert(a.size() == x.size());
  for(int i = 0; i < x.size(); ++i)
    b[i] = a[i] + x[i];

  return;
}

template <typename T, typename S>
void VecDivScalar(std::vector<T> &x, S &alpha)
{
  for(int i = 0; i < x.size(); ++i)
    x[i] = x[i] / alpha;

  return;
}

template <typename S, typename T>
void ScalarVec(S &alpha, std::vector<T> &x)
{
  for(int i = 0; i < x.size(); ++i)
    {
      x[i] = x[i] * alpha;
    }
  return;
}

  
template <typename T>
void Jacobi(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b, const int &MAXITR=10000, const double &tol=10e-6)
{
  std::vector<T> xOld(x.size());
  int itr = 0;
  double error = 100.0;
  for(int itr = 0; itr < MAXITR; ++itr)
    {
      xOld = x;
      for(int i = 0; i < x.size(); ++i)
	{
	  T numerator = 0.0;
	  if(i > 0)
	    numerator += A.lower[i] * xOld[i - 1];
	  if(i >= A.Nx)
	    numerator += A.llower[i] * xOld[i - A.Nx];
	  if(i < A.Nx * A.Nx)
	    numerator += A.upper[i] * xOld[i + 1];
	  if(i < A.Nx * A.Nx - A.Nx)
	    numerator += A.uupper[i] + xOld[i + A.Nx];
	  numerator = -numerator + b[i];
	  x[i] = numerator / A.diag[i];
	}
      
      //if(itr % 10 == 0 && itr > 0)
	error = Norm(x, xOld);

      if(error < tol)
	{
	  std::cerr << "Itr: " << itr << std::endl;
	  return;
	}
    }

  std::cerr << "MAX ITERATIONS EXCEEDED IN JACOBI ITERATIVE (" << MAXITR << ")" << std::endl;
}

template <typename T>
void GaussSeidel(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b, const int &MAXITR=10000, const double &tol=10e-6)
{
  std::vector<T> xOld(x.size());
  double error = 100.0;

  for(int itr = 0; itr < MAXITR; ++itr)
    {
      xOld = x;
      for(int i = 0; i < x.size(); ++i)
	{
	  double numerator = 0.0;
	  if(i > 0)
	    numerator += A.lower[i] * x[i - 1];
	  if(i >= A.Nx)
	    numerator += A.llower[i] * x[i - A.Nx];
	  if(i < A.Nx * A.Nx)
	    numerator += A.upper[i] * xOld[i + 1];
	  if(i < A.Nx * A.Nx - A.Nx)
	    numerator += A.uupper[i] + xOld[i + A.Nx];
	  numerator = -numerator + b[i];
	  x[i] = numerator/A.diag[i];
	}
      
      //if(itr % 10 == 0 && itr > 0)
      if(itr > 0)
	error = Norm(x, xOld);

      if(error < tol)
	{
	  std::cerr << "Itr: " << itr << std::endl;
	  return;
	}
    }
  std::cerr << "MAX ITERATIONS EXCEEDED IN GAUSS-SEIDEL (" << MAXITR << ")" << std::endl; 
}

template <typename T>
void GaussSeidel_SOR(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b, const double &alpha, const int &MAXITR=10000, const double &tol=10e-6)
{
  std::vector<T> xOld(x.size());
  double error = 100.0; 

  for(int itr = 0; itr < MAXITR; ++itr)
    {
      xOld = x;
      for(int i = 0; i < x.size(); ++i)
	{
	  T xOld_i = (1 - alpha) * xOld[i];
	  T numerator = 0.0;
	  
	  if(i > 0)
	    numerator += A.lower[i] * x[i - 1];
	  if(i >= A.Nx)
	    numerator += A.llower[i] * x[i - A.Nx];
	  if(i < A.Nx * A.Nx)
	    numerator += A.upper[i] * xOld[i + 1];
	  if(i < A.Nx * A.Nx - A.Nx)
	    numerator += A.uupper[i] + xOld[i + A.Nx];

	  numerator = alpha * (-numerator + b[i]);
	  x[i] = xOld_i + numerator/A.diag[i];
	}
      //if(itr % 10 == 0 && itr > 0)
	error = Norm(x, xOld);
      
      if(error < tol)
	{
	  std::cerr << "Itr: " << itr << std::endl;
	  return;
	}
    }
  std::cerr << "MAX ITERATIONS EXCEEDED IN GAUSS-SEIDEL SOR (" << MAXITR << ")" << std::endl;
}

template <typename T>
void PCG(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b, const int &MAXITR, const double &tol=10e-6)
{
  std::vector<T> r(x.size(), 0);
  std::vector<T> d(x.size(), 0);
  std::vector<T> q(x.size(), 0);
  std::vector<T> s(x.size(), 0);
  std::vector<T> temp(x.size(), 0);
  Matrix<T> M;
  double del_new = 0.0;
  double del_naught = 0.0;
  double del_old = 0.0;
  double error = 100.0;
  double alpha = 0;
  double beta = 0;
  double e2 = tol*tol;

  // Builds the preconditioner
  M.diag.resize(A.diag.size());
  M.lower.resize(A.lower.size());
  M.llower.resize(A.llower.size());
  M.upper.resize(A.upper.size());
  M.uupper.resize(A.uupper.size());

  for(int i = 0; i < A.diag.size(); ++i)
    {
      M.diag[i] = 1/A.diag[i];
      M.lower[i] = 0;
      M.llower[i] = 0;
      M.upper[i] = 0;
      M.uupper[i] = 0;
    }

  MatVec(A, x, temp);
  VecSubVec(b, temp, r);

  MatVec(M, r, d);
  
  VecVec(r, d, del_new);

  del_naught = del_new;

  for(int itr = 0; itr < MAXITR; ++itr)
    {
      MatVec(A, d, q);
      
      VecVec(d, q, alpha);
      alpha = del_new/alpha;
      
      ScalarVec(alpha, d);
      VecPlusVec(x, d, x);

      if(itr % 50 == 0)
	{
	  MatVec(A, x, temp);
	  VecSubVec(b, temp, r);
	}
      else
	{
	  ScalarVec(alpha, q);
	  VecSubVec(r, q, r);
	}

      MatVec(M, r, s);
      del_old = del_new;
      VecVec(r, s, del_new);
      beta = del_new/del_old;
      ScalarVec(beta, d);
      VecPlusVec(s, d, d);

      if(del_new < e2*del_naught)
	{
	  std::cout << "Itr: " << itr << std::endl;
	  return;
	}
    }
  std::cerr << "MAX ITERATIONS EXCEEDED IN PRECONJUGATE GRADIENT (" << MAXITR << ")" << std::endl;
}

template <typename T>
void GMRES(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b, const int &MAXITR, const double &tol=10e-6)
{
  // www.asc.tuwien.ac.at/~melenk/teach/iterative_SS07/part5.ps
  std::vector<T> r_naught(x.size());
  std::vector<T> v(x.size());
  std::vector<T> w(x.size());
  std::vector<T> temp(x.size());
  double beta = 0.0;

  MatVec(A, x, temp);
  VecSubVec(b, temp, r_naught);

  beta = Norm(r_naught, 0);
  VecDifScalar(v, r_naught);

  std::cerr << "MAX ITERATIONS EXCEEDED IN GMRES (" << MAXITR << ")" << std::endl;
  
}

#endif
