#ifndef QUADTREE_H
#define QUADTREE_H

#include <vector>
#include "planet.h"
#include "vector.h"
#include "objects.h"

template <typename T>
class QuadTree
{
 public:
  QuadTree();
  QuadTree(std::vector<Planet<T> > &bodies, Dimensions &dim);
  ~QuadTree();
  void Insert(Planet<T> &body);
  int maxCount = 4;
 private:
};


template <typename T>
QuadTree<T>::QuadTree()
{
  // Do nothing
}

template <typename T>
QuadTree<T>::QuadTree(std::vector<Planet<T> > &bodies, Dimensions &dim)
{
  for(int i = 0; i < bodies.size(); ++i)
    Insert(bodies[i], i);
  // Do something
}

template <typename T>
void QuadTree<T>::Insert(

#endif
