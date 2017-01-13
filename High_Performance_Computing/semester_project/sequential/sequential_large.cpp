#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include "objects.h"
#include "vector.h"
#include "planet.h"
#include "hermite.h"
//#include "quadtree.h"

using namespace std;

void ReadFile(vector<Planet<double> > &planets, Dimensions &dim);

int main(int argc, char *argv[])
{
  Dimensions dim;

  if(argc == 1)
    dim.numPlanets = 1000;
  else
    dim.numPlanets = atoi(argv[1]);

  double dt = 0.001;
  int years = 2;
  int steps = years/dt;
  steps = 100;//00;

  vector<Planet<double> > planets;

  ReadFile(planets, dim);

  cout << "Start: " << planets[0] << endl;

  Hermite(planets, dt, steps);

  cout << "End: " << planets[0] << endl;

  cout << "Done... Writing Results..." << endl;

  string outfile = to_string(dim.numPlanets);
  outfile.append("done_bodies.dat");

  ofstream output;

  output.open(outfile.c_str());

  for(int i = 0; i < planets.size(); ++i)
    output << planets[i] << endl;

  output.close();
  return 0;
}

void ReadFile(vector<Planet<double> > &planets, Dimensions &dim)
{
  string filename = to_string(dim.numPlanets);
  filename.append("_bodies.dat");
  vector<double> dimensions;
  vector<double> bodies;
  string delimeter = ",";
  string token;
  fstream data;
  size_t pos = 0;
  data.open(filename.c_str());
  string line;

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
