#ifndef _VECTOR_H
#define _VECTOR_H

class Vector
{
 public:
  Vector();
  Vector(double X, double Y, double Z);
  ~Vector();

  Vector operator +(const Vector &rhs);

  void SetX(double X);
  void SetY(double Y);
  void SetZ(double Z);
  double GetX() const;
  double GetY() const;
  double GetZ() const;

 private:
  double x;
  double y;
  double z;
};

#endif
