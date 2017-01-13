#include "vector.h"

Vector::Vector()
{
  x = 0.0;
  y = 0.0;
  z = 0.0;
}

Vector::Vector(double X, double Y, double Z)
{
  x = X;
  y = Y;
  z = Z;
}

Vector::~Vector()
{

}

Vector Vector::operator+(const Vector &rhs)
{
  Vector temp;
  temp.x = GetX() + rhs.GetX();
  temp.y = GetY() + rhs.GetY();
  temp.z = GetZ() + rhs.GetZ();
  return temp;
}

void Vector::SetX(double X)
{
  x = X;
}

void Vector::SetY(double Y)
{
  y = Y;
}

void Vector::SetZ(double Z)
{
  z = Z;
}

double Vector::GetX() const
{
  return x;
}

double Vector::GetY() const
{
  return y;
}

double Vector::GetZ() const
{
  return z;
}
