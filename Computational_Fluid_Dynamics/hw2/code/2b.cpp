#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void WriteFile(vector<double> x, vector<double> u, string time);

int main()
{
  double h = 0.04;
  double dy = 0.0005;
  double dy2 = dy*dy;
  double dt = dy2/20.0;
  int n = h/dy;
  double nu = 0.000217;
  double rho = 800.0;
  double dpdx = 2000.0;
  double beta = 1/rho * dpdx;
  int num_itr = 1.0/dt;
  vector<double> x(n);
  vector<double> u(n);
  double time = 0.0;
  int count = 0;
  int append_count = 0.2/dt;

  u[0] = 40.0;
  x[0] = 0.0;
  for(int i = 1; i < n; ++i)
    {
      x[i] = x[i-1] + dy;
      u[i] = 0.0;
    }

  WriteFile(x, u, to_string(time));

  for(int i = 0; i <= num_itr; ++i)
    {
      time = time + dt;
      for(int j = 1; j < n-1; ++j)
	u[j] = dt * (nu * (u[j+1] - 2*u[j] + u[j-1]) / dy2 - beta) + u[j];
      if(count == append_count)
	{
	  WriteFile(x, u, to_string(time));
	  count = 0;
	}
      ++count;
    }
  return 0;
}

void WriteFile(vector<double> x, vector<double> u, string time)
{
  stringstream ss;
  ss << time[0] << time[1] << time[2] << time[3];
  string true_time = ss.str();
  cout << true_time << " seconds" << endl;
  string outputfile = true_time + "_sec.dat";
  string sys_cmd = "./plot_2.py " + true_time;
  ofstream outfile;
  outfile.open(outputfile.c_str());
  for(int i = 0; i < u.size(); ++i)
    {
      outfile << x[i] << ',' << u[i] << endl;
    }
  outfile.close();

  system(sys_cmd.c_str());

  return;
}
