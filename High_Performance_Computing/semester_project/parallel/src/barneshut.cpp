#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include "mpi.h"
#include "vector.h"
#include "planet.h"
#include "hermite.h"
#include "DendroOctTree.h"

#define HILBERT_ORDERING

using namespace std;

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

void ReadFile(vector<Planet<double> > &planets, Dimensions &dim);

int main(int argc, char *argv[])
{
  int node_id;
  int nodes;
  Dimensions dimension;
  unsigned int dim = 3;
  unsigned int maxDepth = 8;
  unsigned int maxRange = 1<<maxDepth;
  unsigned int maxNumPts = 1;
  double gSize[3] = {1.0, 1.0, 1.0};

  if(argc == 1)
    dimension.numPlanets = 1000;
  else
    dimension.numPlanets = atoi(argv[1]);

  double dt = 0.001;
  int years = 2;
  int steps = years/dt;
  int numPts = 0;
  steps = 100;//00;
  
  vector<Planet<double> > planets;
  vector<double> allPts;
  vector<double> pts;
  vector<ot::TreeNode> tmpNodes;
  vector<ot::TreeNode> linOct;
  vector<ot::TreeNode> balOct;

  MPI::Init(argc, argv);

  nodes = MPI::COMM_WORLD.Get_size();
  node_id = MPI::COMM_WORLD.Get_rank();

  if(node_id == 0)
    {
      ReadFile(planets, dimension);
      numPts = planets.size();
    }

  MPI::COMM_WORLD.Bcast(&numPts, 1, MPI::INT, 0);

  if(numPts % nodes != 0)
    {
      if(node_id == 0)
	cerr << "Number of points isn't divisible by rank. Exiting..." << endl;
      MPI::COMM_WORLD.Abort(-1);
      MPI::Finalize();
      return -1;
    }


  if(node_id == 0)
    {
      double x_scale = dimension.maxX;
      double y_scale = dimension.maxY;
      double z_scale = dimension.maxZ;
      Vector<double> pos;
      for(int i = 0; i < planets.size(); ++i)
	{
	  pos = planets[i].GetPos();
	  double x_s = pos.GetX()/x_scale;
	  double y_s = pos.GetY()/y_scale;
	  double z_s = pos.GetZ()/z_scale;
	  allPts.push_back(x_s);
	  allPts.push_back(y_s);
	  allPts.push_back(z_s);
	}
      numPts = allPts.size()/nodes;
    }

  MPI::COMM_WORLD.Bcast(&numPts, 1, MPI::INT, 0);

  pts.resize(numPts);

  MPI::COMM_WORLD.Scatter(&allPts.front(), numPts/nodes, MPI::DOUBLE, &pts.front(), numPts, MPI::DOUBLE, 0);

  //Converts pts into TreeNode
  for(int i = 0; i < pts.size(); i+=3) {
    if( (pts[i] > 0.0) &&
        (pts[i+1] > 0.0)
        && (pts[i+2] > 0.0) &&
        ( ((unsigned int)(pts[i]*((double)(1u << maxDepth)))) < (1u << maxDepth))  &&
        ( ((unsigned int)(pts[i+1]*((double)(1u << maxDepth)))) < (1u << maxDepth))  &&
        ( ((unsigned int)(pts[i+2]*((double)(1u << maxDepth)))) < (1u << maxDepth)) ) {
      tmpNodes.push_back( ot::TreeNode((unsigned int)(pts[i]*(double)(1u << maxDepth)),
            (unsigned int)(pts[i+1]*(double)(1u << maxDepth)),
            (unsigned int)(pts[i+2]*(double)(1u << maxDepth)),
            maxDepth,dim,maxDepth) );
    }
  }
  pts.clear();
  // Removes the duplicates since points created from distribution
  par::removeDuplicates<ot::TreeNode>(tmpNodes, false, MPI::COMM_WORLD);
  linOct = tmpNodes;
  tmpNodes.clear();
  par::partitionW<ot::TreeNode>(linOct, NULL, MPI::COMM_WORLD);

  pts.resize(3*(linOct.size()));

  for(int i = 0; i < linOct.size(); ++i)
    {
      pts[3*i] = (((double)(linOct[i].getX())) + 0.5)/((double)(1u << maxDepth));
      pts[(3*i)+1] = (((double)(linOct[i].getY())) + 0.5)/((double)(1u << maxDepth));
      pts[(3*i)+2] = (((double)(linOct[i].getZ())) + 0.5)/((double)(1u << maxDepth));
    }

  linOct.clear();

  ot::points2Octree(pts, gSize, linOct, dim, maxDepth, maxNumPts, MPI::COMM_WORLD);
  
  //cout << "Start: " << planets[0] << endl;
  //cout << "hermite" << endl;
  //Hermite(planets, dt, steps);

  //cout << "End: " << planets[0] << endl;

  MPI::Finalize();
  /*
  cout << "Done... Writing Results..." << endl;

  string outfile = to_string(dim.numPlanets);
  outfile.append("done_bodies.dat");

  ofstream output;

  output.open(outfile.c_str());

  for(int i = 0; i < planets.size(); ++i)
    output << planets[i] << endl;

  output.close();
  */
  return 0;
}

void ReadFile(vector<Planet<double> > &planets, Dimensions &dim)
{
  string filename = "data/";
  filename.append(to_string(dim.numPlanets));
  filename.append("_bodies.dat");
  vector<double> dimensions;
  vector<double> bodies;
  string delimeter = ",";
  string token;
  fstream data;
  string line;
  size_t pos = 0;
  bool fileExists = ifstream(filename.c_str()).good();

  if(fileExists == false)
    {
      cerr << "File \"" << filename << "\" does not exist. Exiting..." << endl;
      MPI::COMM_WORLD.Abort(-1);
      MPI::Finalize();
      exit(-1);
    }

  data.open(filename.c_str());
  getline(data, line);
  getline(data, line);
  while((pos = line.find(delimeter)) != string::npos)
    {
      token = line.substr(0, pos);
      dimensions.push_back(stod(token));
      line.erase(0, pos + delimeter.length());
    }

  dim.minX = dimensions[0];
  dim.maxX = dimensions[1];
  dim.minY = dimensions[2];
  dim.maxY = dimensions[3];
  dim.minZ = dimensions[4];
  dim.maxZ = stod(line);

  while(getline(data, line))
    {
      pos = 0;
      bodies.clear();
      while((pos = line.find(delimeter)) != string::npos)
	{
	  token = line.substr(0, pos);
	  bodies.push_back(stod(token));
	  line.erase(0, pos + delimeter.length());
	}
      bodies.push_back(stod(line));
      Planet<double> temp(bodies[1], bodies[2], bodies[3], bodies[0]);
      planets.push_back(temp);
    }

  data.close();

  cout << "Finished Reading Data" << endl;

  return;
}
