#ifndef DIMENSIONS_H
#define DIMENSIONS_H
#include <vector>

template <typename T>
struct Dimensions
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
  std::vector<T> x;
  std::vector<T> y;
};

template <typename T>
struct Flow
{
  std::vector<T> psi;
  std::vector<T> omega;
  std::vector<T> u;
  std::vector<T> v;
};

#endif
