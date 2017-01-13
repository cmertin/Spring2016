#ifndef MULTIGRID_H
#define MULTIGRID_H
#include <vector>
#include <cmath>
#include "Matrix.h"

struct MultiGrid
{
  double omega;
  const int MAXLEVEL;
  int nu_1;
  int nu_2;
};

template<typename T>
void RichardsonSmoothing(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b, const double &omega)
{
  std::vector<T> r(x.size());
  MatVec(A, x, r);
  VecSubVec(b, r, r);
  for(int i = 0; i < b.size(); ++i)
    b[i] = b[i] + omega * r[i];

  return;
}

template<typename T>
void JacobiSmoothing(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b, const double &omega)
{
  std::vector<T> r(x.size());
  std::vector<T> b_tmp(x.size());
  Matrix<T> M;

  M.diag.resize(A.diag.size());
  M.upper.resize(A.upper.size());
  M.lower.resize(A.lower.size());
  M.uupper.resize(A.uupper.size());
  M.llower.resize(A.llower.size());
  
  for(int i = 0; i < x.size(); ++i)
    {
      M.diag[i] = 1/A.diag[i];
      M.upper[i] = 0;
      M.lower[i] = 0;
      M.uupper[i] = 0;
      M.llower[i] = 0;
    }
  
  MatVec(A, x, r);
  VecSubVec(b, r, r);

  Jacobi(M, b_tmp, r);

  for(int i = 0; i < b.size(); ++i)
    b[i] = x[i] + omega * b_tmp;

  return;
}

template <typename T>
void SmoothGaussSeidel(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b)
{
  std::vector<T> r(x.size());
  std::vector<T> b_tmp(x.size());
  Matrix<T> M;

  M.diag.resize(A.diag.size());
  M.upper.resize(A.upper.size());
  M.lower.resize(A.lower.size());
  M.uupper.resize(A.uupper.size());
  M.llower.resize(A.llower.size());
  
  for(int i = 0; i < x.size(); ++i)
    {
      M.diag[i] = A.diag[i];
      M.upper[i] = 0;
      M.lower[i] = A.lower[i];
      M.uupper[i] = 0;
      M.llower[i] = A.llower[i];
    }

  GaussSeidel(M, b_tmp, r);
  for(int i = 0; i < b.size(); ++i)
    b[i] = x[i] + b_tmp[i];

  return;
}

/*
  r = b - A*x
  D = diag(diag(A))
  E = D - tril(A) 
  M = (D - omega*E)/omega
  y = x + M\r
 */
template <typename T>
void SmoothSOR(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b, const double &omega)
{
  
}

template <typename T>
void RestrictionMatrix(std::vector<std::vector<T> > &R, const int &N)
{
  R.resize(N);
  int nr = (sqrt(N)-1)/2;
  for(int i = 0; i < R.size(); ++i)
    R[i].resize(nr);

  for(int i = 0; i < nr; ++i)
    {
      for(int j = 0; j < nr; ++j)
	{
	  R[i*nr+j][(2*i+1)*(2*nr+1)+2*(j)-1] = 4.0/16;
	  R[i*nr+j][(2*i+1)*(2*nr+1)+2*(j)] = 2.0/16;
	  R[i*nr+j][(2*i+1)*(2*nr+1)+2*(j)+1] = 2.0/16;
	  R[i*nr+j][(2*i+1)*(2*nr+1)+2*(j)+(2*nr+1)-1] = 2.0/16;
	  R[i*nr+j][(2*i+1)*(2*nr+1)+2*(j)-(2*nr+1)-1] = 2.0/16;
	  R[i*nr+j][(2*i+1)*(2*nr+1)+2*j+1+(2*nr+1)-1] = 1.0/16;
	  R[i*nr+j][(2*i+1)*(2*nr+1)+2*j+1-(2*nr+1)-1] = 1.0/16;
	  R[i*nr+j][(2*i+1)*(2*nr+1)+2*j+(2*nr+1)-1-1] = 1.0/16;
	  R[i*nr+j][(2*i+1)*(2*nr+1)+2*j-(2*nr+1)-1-1] = 1.0/16;
	}
    }
}

template <typename T, template S>
void InterpolationMatrix(std::vector<std::vector<T> > &R, std::vector<std::vector<T> > &P, const S &alpha)
{
  int oldRow = P.size();
  int oldCol = P[0].size();

  P.resize(oldCol);
  for(int i = 0; i < P.size(); ++i)
    P[i].resize(oldRow);
  
  for(int i = 0; i < R.size(); ++i)
    {
      for(int j = 0; j < R[i].size(); ++j)
	swap(R[i][j],P[j][i]);
    }

  for(int i = 0; i < P.size(); ++i)
    {
      for(int j = 0; j < P[i].size(); ++j)
	P[i][j] = alpha * P[i][j];
    }

  return;
}

template <typename T>
void MultiGrid(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b, const int &level, MultiGrid &settings)
{
  std::vector<T> rh(x.size());
  std::vector<T> rH(x.size());
  std::vector<T> temp(x.size());
  std::vector<std::vector<T> > IhH;
  std::vector<std::vector<T> > IHh;
  
  
  for(int i = 0; i < nu_1; ++i)
    SmoothGaussSeidel(A, x, b);

  MatVec(A, x, temp);
  VecSubVec(b, temp, rh);

  RestrictionMatrix(R, A.diag.size());
  InterpolationMatrix(R, P, 4);

  if(level == 0)
    deltaH = 
  
}

template <typename T>
void MultiGrid(Matrix<T> &A, std::vector<T> &x, std::vector<T> &b, const int &MAXITR, MultiGrid &settings, const double &tol=10e-6)
{
  MGM(A, x, b, l, MultiGrid &settings);
}

#endif
