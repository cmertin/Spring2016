#ifndef DERIVATIVE_H
#define DERIVATIVE_H

template <typename T>
T inline CentralDifference(T &x, T &h, T (*fn)(T X))
{
  return (fn(x+h) - fn(x-h))/(2*h);
}

// Returns the approximate second order central difference
template <typename T>
T inline CentralDifference_2(T &x, T &h, T (*fn)(T X))
{
  T numerator = fn(x + h) - 2 * fn(x) + fn(x - h);
  T denominator = h*h;
  return numerator/denominator;
}

// Returns the approximate first order forward difference
template <typename T>
T inline ForwardDifference(T &x, T &h, T (*fn)(T X))
{
  return (fn(x+h) - fn(x))/h;
}

// Returns the approximate second order forward difference
template <typename T>
T inline ForwardDifference_2(T &x, T &h, T (*fn)(T X))
{
  T numerator = fn(x + 2*h) - 2*fn(x + h) + fn(x);
  T denominator = h*h;
  return numerator/denominator;
}

#endif
