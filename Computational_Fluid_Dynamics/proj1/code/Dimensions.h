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
  double U_left;
  double U_right;
  double U_bottom;
  double U_top;
  std::vector<std::vector<T> > P;
  std::vector<std::vector<T> > u;
  std::vector<std::vector<T> > ut;
  std::vector<std::vector<T> > v;
  std::vector<std::vector<T> > vt;
  std::vector<std::vector<T> > omega;
};

#endif
