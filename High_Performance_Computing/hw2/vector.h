#ifndef _VECTOR_H
#define _VECTOR_H

#include <iostream>
#include <cmath>

template <typename T>
class Vector
{
 public:
  Vector();
  Vector(T x, T y, T z);
  ~Vector();

  friend std::ostream &operator<<(std::ostream &out, const Vector<T> &vec)
  {
    out << '(' << vec.x << ", " << vec.y << ", " << vec.z << ')';
    return out;
  }
  // Vector assignment/equate
  Vector<T> &operator=(const Vector<T> &rhs);
  // Vector operations
  Vector<T> operator+(const Vector<T> &rhs) const;
  Vector<T> &operator+=(const Vector<T> &rhs);
  Vector<T> operator-(const Vector<T> &rhs) const;
  Vector<T> operator*(const Vector<T> &rhs) const;
  Vector<T> operator/(const Vector<T> &rhs) const;

  T Magnitude();

  void SetX(T x);
  void SetY(T y);
  void SetZ(T z);
  T GetX() const;
  T GetY() const;
  T GetZ() const;

 private:
  T x;
  T y;
  T z;
};

template <typename T>
Vector<T>::Vector()
{
  this->x = 0.0;
  this->y = 0.0;
  this->z = 0.0;
}

template <typename T>
Vector<T>::Vector(T x, T y, T z)
{
  this->x = x;
  this->y = y;
  this->z = z;
}

template <typename T>
Vector<T>::~Vector()
{
  // do nothing
}

template <typename T>
Vector<T> &Vector<T>::operator=(const Vector<T> &rhs)
{
  this->x = rhs.GetX();
  this->y = rhs.GetY();
  this->z = rhs.GetZ();
  return *this;
}

template <typename T>
Vector<T> Vector<T>::operator+(const Vector<T> &rhs) const
{
  Vector<T> result;
  result.x = this->x + rhs.GetX();
  result.y = this->y + rhs.GetY();
  result.z = this->z + rhs.GetZ();
  return result;
}

template <typename T>
Vector<T> &Vector<T>::operator+=(const Vector<T> &rhs)
{
  this->x = this->x + this->x;
  this->y = this->y + this->y;
  this->z = this->z + this->z;
  return *this;
}

template <typename T>
Vector<T> Vector<T>::operator-(const Vector<T> &rhs) const
{
  Vector<T> result;
  result.x = this->x - rhs.GetX();
  result.y = this->y - rhs.GetY();
  result.z = this->z - rhs.GetZ();

  return result;
}

template <typename T>
Vector<T> Vector<T>::operator*(const Vector<T> &rhs) const
{
  Vector<T> result;
  result.x = this->x * rhs.GetX();
  result.y = this->y * rhs.GetY();
  result.z = this->z * rhs.GetZ();

  return result;
}

template <typename T>
Vector<T> Vector<T>::operator/(const Vector<T> &rhs) const
{
  Vector<T> result;
  result.x = this->x / rhs.GetX();
  result.y = this->y / rhs.GetY();
  result.z = this->z / rhs.GetZ();

  return result;
}

template <typename T>
T Vector<T>::Magnitude()
{
  return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}

template <typename T>
void Vector<T>::SetX(T x)
{
  this->x = x;
}

template <typename T>
void Vector<T>::SetY(T y)
{
  this->y = y;
}

template <typename T>
void Vector<T>::SetZ(T z)
{
  this->z = z;
}

template <typename T>
T Vector<T>::GetX() const
{
  return this->x;
}

template <typename T>
T Vector<T>::GetY() const
{
  return this->y;
}

template <typename T>
T Vector<T>::GetZ() const
{
  return this->z;
}

#endif
