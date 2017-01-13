#ifndef _POINTS_H
#define _POINTS_H

#include <cmath>
#include <iostream>

template <typename T>
class Point
{
 public:
  Point(); // Default Constructor
  Point(T x, T y); // Constructor
  ~Point(); // Destructor
  // Copy Constructor

  // Comparison Operators
  bool operator <(const Point<T> &rhs) const;
  bool operator >(const Point<T> &rhs) const;
  friend std::ostream& operator<<(std::ostream &out, const Point<T> &point)
  {
    out << '(' << point.x << ", " << point.y << ')';
    return out;
  }

  // Normal Function Calls
  T GetX();
  T GetY();
  T Magnitude();
  void SetX(T X);
  void SetY(T Y);
  void Print();
  
  // Standard
  bool Equal(T a, T b, float tol=0.00004) const;
  
 private:
  T x;
  T y;
};

template <typename T>
Point<T>::Point()
{
  this->x = 0.0;
  this->y = 0.0;
}

template <typename T>
Point<T>::Point(T x, T y)
{
  this->x = x;
  this->y = y;
}

template <typename T>
Point<T>::~Point()
{
  // Nothing to do
}

template <typename T>
bool Point<T>::operator <(const Point &rhs) const
{
  if(this->y > rhs.y && !Equal(this->y, rhs.y))
    return false;
  else
    {
      if(Equal(this->y, rhs.y) && this->x > rhs.x)
	return false;
      else
	return true;
    }
}

template <typename T>
bool Point<T>::operator >(const Point<T> &rhs) const
{
  if(this->y > rhs.y && !Equal(this->y, rhs.y))
    return true;
  else
    {
      if(Equal(this->y, rhs.y) && this->x > rhs.x)
	return true;
      else
	return false;
    }
}


template <typename T>
T Point<T>::GetX()
{
  return this->x;
}

template <typename T>
T Point<T>::GetY()
{
  return this->y;
}

template <typename T>
T Point<T>::Magnitude()
{
  return sqrt(this->x*this->x + this->y*this->y);
}

template <typename T>
void Point<T>::SetX(T X)
{
  this->x = X;
  return;
}

template <typename T>
void Point<T>::SetY(T X)
{
  this->y = X;
  return;
}

template <typename T>
void Point<T>::Print()
{
  std::cout << '(' << this->x << ", " << this->y << ')' << std::endl;
}

template <typename T>
bool Point<T>::Equal(T a, T b, float tol) const
{
  if(std::abs(a - b) < tol)
    return true;
  else
    return false;
}

#endif
