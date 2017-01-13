#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

struct Geometry
{
  double minX;
  double maxX;
  double minY;
  double maxY;
  double dx;
  double dy;
  double dt;
  double U;
  double Re;
  double nu;
  int Nx;
  int Ny;
  vector<double> x;
  vector<double> y;
};

bool Error(vector<double> &omega, vector<double> &omegaOld, Geometry &geom, double &tol);
void UpdateVelocities(vector<double> &psi, vector<double> &u, vector<double> &v);
inline double Omega_fn(vector<double> &omega, vector<double> &psi, Geometry &geom, int &i, int &j);
inline double Psi_fn(vector<double> &omega, vector<double> &psi, Geometry &geom, int &i, int &j);

int main()
{
  // Initialize variables
  Geometry geom;
  geom.minX = 0;
  geom.maxX = 1.0;
  geom.minY = 0;
  geom.maxY = 1.0;
  geom.Nx = 50;
  geom.Ny = 50;
  geom.Re = 100;
  geom.nu = 0.1;
  geom.dt = 10e-4;
  geom.dx = (geom.maxX - geom.minX)/(geom.Nx - 1);
  geom.dy = (geom.maxY - geom.minY)/(geom.Ny - 1);
  geom.U = (geom.Re * geom.nu)/(geom.maxX - geom.minX);

  vector<double> omega(geom.Nx * geom.Ny, 0);
  vector<double> omegaOld(geom.Nx * geom.Ny, 0);
  vector<double> psi(geom.Nx * geom.Ny, 0);
  vector<double> psiOld(geom.Nx * geom.Ny, 0);
  vector<double> u(geom.Nx * geom.Ny, 0);
  vector<double> v(geom.Nx * geom.Ny, 0);
  geom.x.resize(geom.Nx, 0);
  geom.y.resize(geom.Nx, 0);
  double tol = 10e-6;
  int count = 0;
  bool error = false;

  // Initialize parameters
  for(int i = 1; i < geom.Nx; ++i)
    geom.x[i] = geom.x[i-1] + geom.dx;

  for(int i = 1; i < geom.Ny; ++i)
    geom.y[i] = geom.y[i-1] + geom.dy;

  for(int i = 0; i < geom.Nx; ++i)
    u[i + (geom.Ny-1)*geom.Nx] = geom.U;
  
  for(int j = 1; j < geom.Ny-1; ++j)
    {
      double dx_ = 2.0*geom.dx;
      double dy_ = 2.0*geom.dy;
      double dv_dx;
      double du_dy;
      for(int i = 1; i < geom.Nx-1; ++i)
	{
	  int index = i + j*geom.Nx;
	  dv_dx = (v[i+1 + j*geom.Nx] - v[i-1 + j*geom.Nx])/dx_;
	  du_dy = (u[i + (j+1)*geom.Nx] - u[i + (j-1)*geom.Nx])/dy_;
	  omega[index] = dv_dx - du_dy;
	}
    }

  // Perform the calculations
  do
    {
      omegaOld = omega;
      psiOld = psi;
      double dx2 = geom.dx * geom.dx;
      for(int j = 1; j < geom.Ny-1; ++j)
	{
	  for(int i = 1; i < geom.Nx-1; ++i)
	    {
	      int index = i + j*geom.Nx;
	      psi[index] = Psi_fn(omegaOld, psiOld, geom, i, j);
	      omega[index] = Omega_fn(omegaOld, psiOld, geom, i, j);
	      psiOld[index] = psi[index];
	      omegaOld[index] = omega[index];
	      omega[0 + j*geom.Nx] = -2 * (psi[1 + j * geom.Nx])/dx2;
	      omega[i + 0*geom.Nx] = -2 * (psi[i * 1 * geom.Nx])/dx2;
	      omega[geom.Nx-1 + j * geom.Nx] = -2 * (psi[geom.Nx-2 + j*geom.Nx])/dx2;
	      omega[i + (geom.Ny-1)*geom.Nx] = -2 * (psi[i + (geom.Ny - 1)*geom.Nx] + geom.U * geom.dx)/dx2;
	    }
	}
      ++count;
      error = Error(omega, psi, geom, tol);
    }while(error != true);

  cout << "Iterations: " << count << endl;

  for(int j = 0; j < geom.Ny; ++j)
    {
      double dx_ = 2*geom.dx;
      double dy_ = 2*geom.dy;
      int N = geom.Nx;
      for(int i = 0; i < geom.Nx; ++i)
	{
	  int index = i + j*N;
	  double dpsi_dy = (psi[i + (j+1)*N] - psi[i + (j-1)*N])/dx_;
	  double dpsi_dx = (psi[i+1 + j*N] - psi[i-1 + j*N])/dy_;
	  if(j == geom.Nx - 1) // Top boundary
	    dpsi_dy = geom.U;
	  else if(j == 0) // Bottom boundary
	    dpsi_dy = 0;
	  else if(i == 0) // Left boundary
	    dpsi_dx = 0;
	  else if(i == geom.Nx - 1)
	    dpsi_dx = 0;
	  u[index] = dpsi_dy;
	  v[index] = -dpsi_dx;
	}
    }

  for(int j = 0; j < geom.y.size(); ++j)
    {
      for(int i = 0; i < geom.y.size(); ++i)
	{
	  int index = i + j*geom.Nx;
	  cout << geom.x[i] << ',' << geom.y[i] << ',' << omega[index] << ',' << psi[index] << ',' << u[index] << ',' << v[index] << endl;
	}
    }
  return 0;
}

void UpdateVelocities(vector<double> &psi, vector<double> &u, vector<double> &v)
{
  
}

bool Error(vector<double> &omega, vector<double> &psi, Geometry &geom, double &tol)
{
  vector<double> diff(omega.size());

  for(int j = 1; j < geom.Ny-1; ++j)
    {
      for(int i = 1; i < geom.Nx-1; ++i)
	{
	  int index = i + j * geom.Nx;
	  diff[index] = abs(omega[index] - Omega_fn(omega, psi, geom, i, j));
	}
    }

  double max = *max_element(diff.begin(), diff.end());

  if(max < tol)
    return true;
  else
    return false;
}

inline double Omega_fn(vector<double> &omega, vector<double> &psi, Geometry &geom, int &i, int &j)
{
  int N = geom.Nx;
  double firstTerm = omega[i+1 + j*N] + omega[i-1 + j*N] + omega[i + (j+1)*N] + omega[i + (j-1)*N];
  double secTerm = (psi[i + (j+1)*N] - psi[i + (j-1)*N])*(omega[i+1 + j*N] - omega[i-1 + j*N]);
  double thirTerm = (psi[i+1 + j*N] - psi[i-1 + j*N])*(omega[i + (j+1)*N] - omega[i + (j-1)*N]);
  secTerm = (geom.Re/4.0) * (secTerm - thirTerm);
  return (1.0/4.0) * (firstTerm - secTerm);
}

inline double Psi_fn(vector<double> &omega, vector<double> &psi, Geometry &geom, int &i, int &j)
{
  int N = geom.Nx;
  double dx2 = geom.dx * geom.dx;
  return 0.25 * (psi[i+1 + j*N] + psi[i-1 + j*N] + psi[i + (j+1)*N] + psi[i + (j-1)*N] + dx2 * omega[i + j*N]);
}
