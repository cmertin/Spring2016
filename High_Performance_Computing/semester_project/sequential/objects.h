#ifndef OBJECTS_H
#define OBJECTS_H

struct Dimensions
{
  double minX;
  double maxX;
  double minY;
  double maxY;
  double minZ;
  double maxZ;
  int numPlanets;
};

struct NodeDim
{
  double x;
  double y;
};

class Node
{
 public:
  double x,y;
  double radius;
  int id; // Unique ID of node
  int inNode = 0;
  Node *nwChild = NULL;
  Node *neChild = NULL;
  Node *swChild = NULL;
  Node *seChild = NULL;
  Node *parent = NULL;
  /*
    for octtree:
    front, back, left, right, top, bottom
   */
  Node(NodeDim &dim, int &id)
    {
      x = dim.x;
      y = dim.y;
      this->id = id;
    }
};

#endif
