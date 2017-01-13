#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include "planet.h"

using namespace std;

void WriteBodies(vector<Planet<double> > &planets, double minX, double maxX);

int main(int argc, char *argv[])
{
  int maxNum = 10000;
  // Random number generator
  double minX = 0;//-10000
  double maxX = 10;//10000
  double minMass = 10e-8;
  double maxMass = 10;
  random_device rd;
  default_random_engine generator(rd()); // rd() provides a random seed
  uniform_real_distribution<double> posGen(minX, maxX);
  uniform_real_distribution<double> massGen(minMass, maxMass);

  vector<Planet<double> > planets;

  for(int i = 1; i <= maxNum; ++i)
    {
      double mass = massGen(generator);
      double x = posGen(generator);
      double y = posGen(generator);
      double z = posGen(generator);
      Planet<double> temp(x, y, z, mass);
      planets.push_back(temp);
      
      if(i == 10)
	WriteBodies(planets, minX, maxX);
      if(i == 100)
	WriteBodies(planets, minX, maxX);
      
      if(i % 1000 == 0)
	WriteBodies(planets, minX, maxX);
    }

  
  return 0;
}

void WriteBodies(vector<Planet<double> > &planets, double minX, double maxX)
{
  string filename = "data/";
  filename.append(to_string(planets.size()));
  filename.append("_bodies.dat");
  ofstream output;
  output.open(filename.c_str());
  output << planets.size() << endl;
  output << minX << ',' << maxX << ',' << minX << ',' << maxX << ',' << minX << ',' << maxX << endl;

  for(int i = 0; i < planets.size(); ++i)
    output << planets[i] << endl;
  output.close();
}
    
