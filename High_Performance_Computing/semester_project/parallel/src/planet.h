#ifndef PLANET_H
#define PLANET_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include "vector.h"

template <typename T>
class Planet
{
 public:
  static constexpr const T G = 6.67408e-11;
  Planet();
  Planet(T x, T y, T z, T mass);
  ~Planet();

  friend std::ostream &operator<<(std::ostream &out, const Planet<T> &planet)
  {
    out << std::setprecision(5) << std::scientific << planet.mass << ',' << std::fixed 
	<< planet.pos.GetX() << ',' << planet.pos.GetY() << ',' << planet.pos.GetZ();
    return out;
  }
  Vector<T> Force(Planet<T> &rhs);

  Vector<T> GetPos();
  Vector<T> GetVel();
  Vector<T> GetAcc();
  Vector<T> GetJerk();
  Vector<T> GetPosOld();
  Vector<T> GetVelOld();
  Vector<T> GetAccOld();
  Vector<T> GetJerkOld();
  T GetMass();

  void SetPos(Vector<T> &rhs);
  void SetPosOld(Vector<T> &rhs);
  void SetAcc(Vector<T> &rhs);
  void SetAccOld(Vector<T> &rhs);
  void SetVel(Vector<T> &rhs);
  void SetVelOld(Vector<T> &rhs);
  void SetJerk(Vector<T> &rhs);
  void SetJerkOld(Vector<T> &rhs);
  void SetMass(T mass);

 private:
  Vector<T> pos;
  Vector<T> pos_old;
  Vector<T> vel;
  Vector<T> vel_old;
  Vector<T> acc;
  Vector<T> acc_old;
  Vector<T> jerk;
  Vector<T> jerk_old;
  T mass;
};

template <typename T>
Planet<T>::Planet()
{
  Vector<T> zero;
  this->pos = zero;
  this->pos_old = zero;
  this->vel = zero;
  this->vel_old = zero;
  this->acc = zero;
  this->acc_old = zero;
  this->jerk = zero;
  this->jerk_old = zero;
  this->mass = 0;
}

template <typename T>
Planet<T>::Planet(T x, T y, T z, T mass)
{
  Vector<T> A(x, y, z);
  Vector<T> zero;
  this->pos = A;
  this->pos_old = zero;
  this->vel = zero;
  this->vel_old = zero;
  this->acc = zero;
  this->acc_old = zero;
  this->jerk = zero;
  this->jerk_old = zero;
  this->mass = mass;
}

template <typename T>
Planet<T>::~Planet()
{

}

template <typename T>
Vector<T> Planet<T>::Force(Planet<T> &rhs)
{
  Vector<T> vec = this->pos - rhs.pos;
  double coeff = - G * GetMass() * rhs.GetMass()/pow(vec.Magnitude(),3);
  return vec * coeff;
}

template <typename T>
Vector<T> Planet<T>::GetPos()
{
  return this->pos;
}

template <typename T>
Vector<T> Planet<T>::GetPosOld()
{
  return this->pos_old;
}

template <typename T>
Vector<T> Planet<T>::GetVel()
{
  return this->vel;
}

template <typename T>
Vector<T> Planet<T>::GetVelOld()
{
  return this->vel_old;
}

template <typename T>
Vector<T> Planet<T>::GetAcc()
{
  return this->acc;
}

template <typename T>
Vector<T> Planet<T>::GetAccOld()
{
  return this->acc_old;
}

template <typename T>
Vector<T> Planet<T>::GetJerk()
{
  return this->jerk;
}

template <typename T>
Vector<T> Planet<T>::GetJerkOld()
{
  return this->jerk_old;
}

template <typename T>
T Planet<T>::GetMass()
{
  return this->mass;
}

template <typename T>
void Planet<T>::SetPos(Vector<T> &rhs)
{
  this->pos = rhs;
}

template <typename T>
void Planet<T>::SetPosOld(Vector<T> &rhs)
{
  this->pos_old = rhs;
}

template <typename T>
void Planet<T>::SetVel(Vector<T> &rhs)
{
  this->vel = rhs;
}

template <typename T>
void Planet<T>::SetVelOld(Vector<T> &rhs)
{
  this->vel_old = rhs;
}

template <typename T>
void Planet<T>::SetAcc(Vector<T> &rhs)
{
  this->acc = rhs;
}

template <typename T>
void Planet<T>::SetAccOld(Vector<T> &rhs)
{
  this->acc_old = rhs;
}

template <typename T>
void Planet<T>::SetJerk(Vector<T> &rhs)
{
  this->jerk = rhs;
}

template <typename T>
void Planet<T>::SetJerkOld(Vector<T> &rhs)
{
  this->jerk_old = rhs;
}

template <typename T>
void Planet<T>::SetMass(T mass)
{
  this->mass = mass;
}

#endif
