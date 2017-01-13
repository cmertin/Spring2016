#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include "vector.h"
#include "planet.h"

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
  Dimensions dim;

  if(argc == 1)
    dim.numPlanets = 1000;
  else
    dim.numPlanets = atoi(argv[1]);
  
  static constexpr const double G = 4 * M_PI * M_PI;
  Vector<double> zero;
  double dt = 0.001;
  int years = 2;
  int steps = years/dt;
  steps = 100;//00;

  vector<Planet<double> > planets;//(dim.numPlanets);

  ReadFile(planets, dim);

  /*
  double e_dist = 1.0;//149.6*10e6;
  double m_dist = 1.38;//227.9*10e6;
  double e23 = sqrt(e_dist * e_dist / 2.0);
  double m23 = -sqrt(m_dist * m_dist / 2.0);
  vector<Planet<double> > planets(3);
  Planet<double> Sun(0, 0, 0, 1.0);
  Planet<double> Earth(e23, e23, 0, 3.003e-6);
  Planet<double> Mars(m23, m23, 0, 3.213e-7);
  vector<Planet<double> > pos;

  Vector<double> test(0.4, -.4, 0);
  Vector<double> test2(0.5, 1.3, 3.7);
  */
  Vector<double> r;
  Vector<double> r3;
  Vector<double> r5;
  Vector<double> v;
  Vector<double> a;
  Vector<double> ai;
  Vector<double> aj;
  Vector<double> jk;
  Vector<double> jki;
  Vector<double> jkj;
      
  /*
  // Initialize some things
  for(int i = 0; i < planets.size(); ++i)
    {
      double phi = i * 2 * M_PI/3.0;
      double vabs = 1.0/sqrt(sqrt(3));
      Vector<double> vel(-vabs * sin(phi) + 0.001, vabs*cos(phi), -vabs*cos(phi) * sin(phi));
      planets[i].SetVel(vel);
    }
  */
  cout << "Start: " << planets[0] << endl;
  
  //cout << planets[2].GetVel() << endl;

  // Do the time steps
  for(int itr = 0; itr < steps; ++itr)
    {
      for(int i = 0; i < planets.size(); ++i)
	{
	  r = planets[i].GetPos();
	  v = planets[i].GetVel();
	  a = planets[i].GetAcc();
	  jk = planets[i].GetJerk();
	  planets[i].SetPosOld(r);
	  planets[i].SetVelOld(v);
	  planets[i].SetAccOld(a);
	  planets[i].SetJerkOld(jk);
	  v += a * dt + jk * (dt*dt)/2.0;
	  r += v * dt + a * ((dt*dt)/2.0) + jk*((dt*dt*dt)/6.0);
	  //if(i == 1)
	    //cout << r << endl;
	  //cout << r << '\t' << planets[i].GetPosOld() << endl;
	  planets[i].SetPos(r);
	  planets[i].SetVel(v);
	}
      for(int i = 0; i < planets.size(); ++i)
	{
	  a = zero;
	  jk = zero;
	  for(int j = 0; j < planets.size(); ++j)
	    {
	      if(i == j)
		continue;
	      Vector<double> rji = planets[j].GetPos() - planets[i].GetPos();
	      Vector<double> vji = planets[j].GetVel() - planets[i].GetVel();
	      double rji3 = pow(rji.Magnitude(),3);
	      double rji5 = pow(rji.Magnitude(),5);
	      double rv = rji * vji;
	      a += rji * planets[j].GetMass()/rji3;
	      //cout << a << endl;
	      //cout << rji.Magnitude() << endl;
	      jk += (vji/rji3 -  rji * ((3/rji5) * rv)) * planets[j].GetMass();
	    }
	  a = a * G;
	  jk = jk * G;
	  //if(i == 1)
	    //cout << a << endl;
	  planets[i].SetAcc(a);
	  planets[i].SetJerk(jk);
	}

      for(int i = 0; i < planets.size(); ++i)
	{
	  v = planets[i].GetVelOld() + (planets[i].GetAcc() + planets[i].GetAccOld()) * (dt/2.0)
	    + (planets[i].GetJerk() - planets[i].GetJerkOld()) * (dt * dt) / 12.0;
	  r = planets[i].GetPos() + (planets[i].GetVel() + planets[i].GetVelOld()) * (dt/2.0)
	    + (planets[i].GetAcc() - planets[i].GetAccOld()) * (dt * dt) / 12.0;
	  planets[i].SetVel(v);
	  planets[i].SetPos(r);
	  //if(i == 1)
	  //cout << v << endl;
	}
    }
  cout << "End: " << planets[0] << endl;

  cout << "Done... Writing Results..." << endl;

  string outfile = to_string(dim.numPlanets);
  outfile.append("done_bodies.dat");

  ofstream output;

  output.open(outfile.c_str());

  for(int i = 0; i < planets.size(); ++i)
    output << planets[i] << endl;

  output.close();
  /*
  string sunfile = "sun.dat";
  string earthfile = "earth.dat";
  string marsfile = "mars.dat";

  ofstream sunout;
  ofstream earthout;
  ofstream marsout;

  sunout.open(sunfile.c_str());
  earthout.open(earthfile.c_str());
  marsout.open(marsfile.c_str());

  for(int i = 0; i < pos.size(); i += 3)
    {
      sunout << pos[i] << endl;
      earthout << pos[i+1] << endl;
      marsout << pos[i+2] << endl;
    }

  sunout.close();
  earthout.close();
  marsout.close();
  */
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
