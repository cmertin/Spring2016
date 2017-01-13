#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <algorithm>
#include "derivative.h"

using namespace std;

template <typename T>
struct errors
{
  vector<T> xVals;
  vector<T> derExact;
  vector<T> derApprox;
  vector<T> errorNorms;
  vector<T> errorsX;
  T h;
  string derLabel; // CD, FWD, CD2
  string derType; // fnP, fnPP
};

double fn(double x);
double fn_prime(double x);
double fn_pprime(double x);
double RelError(double exact, double approx);
double Error(double exact, double approx);
double ErrorNorm(vector<double> exact, vector<double> approx);
template <typename T>
void WriteFile(vector<T> x, vector<T> y, T h, string filename="", bool exact=false);
template <typename T>
void WriteErrors(errors<T> err);

int main(int argc, char *argv[])
{
  double x = 0.5;
  vector<double> h = {0.1, 0.01, 0.001};
  vector<string> labels = {"CD", "FWD", "CD2"};
  vector<string> exactLabels = {"fn_exact", "fnP_exact", "fnPP_exact"};
  vector<string> errorLabels = {"fnP", "fnPP"};
  vector<double> xVals;
  vector<double> derApprox;
  vector<double> fnReal;
  vector<double> derExact;
  double xMin = 0.0;
  double xMax = 1.0;

  // Approximating the derivatives
  for(int lbl = 0; lbl < labels.size(); ++lbl)
    {
      for(int i = 0; i < h.size(); ++i)
	{
	  xVals.clear();
	  derApprox.clear();
	  double steps = (xMax - xMin)/h[i];
	  for(int j = 0; j <= steps; ++j)
	    {
	      double x = xMin + h[i] * j;
	      double fn_p = 0;
	      if(lbl == 0) // Central Difference fn'
		fn_p = CentralDifference(x, h[i], fn);
	      else if(lbl == 1) // Forward Difference fn'
		fn_p = ForwardDifference(x, h[i], fn);
	      else if(lbl == 2) // Central Difference fn''
		fn_p = CentralDifference_2(x, h[i], fn);
	      xVals.push_back(x);
	      derApprox.push_back(fn_p);
	    }
	  WriteFile(xVals, derApprox, h[i], labels[lbl]);
	}
    }

  // Exact values of function and derivatives
  for(int i = 0; i < exactLabels.size(); ++i)
    {
      double stepSize = 0.0001;
      double steps = (xMax - xMin)/stepSize;
      xVals.clear();
      fnReal.clear();
      for(int j = 0; j <= steps; ++j)
	{
	  double x = xMin + stepSize * j;
	  double y = 0;
	  if(i == 0) // Original Function
	    y = fn(x);
	  else if(i == 1)
	    y = fn_prime(x);
	  else if(i == 2)
	    y = fn_pprime(x);
	  xVals.push_back(x);
	  fnReal.push_back(y);
	}
      WriteFile(xVals, fnReal, stepSize, exactLabels[i], true);
    }

  errors<double> err;

  // Error of approximation
  // fn'
  for(int i = 0; i < labels.size()-1; ++i)
    {
      err.derLabel = labels[i];
      err.derType = "fnP";
      for(int j = 0; j < h.size(); ++j)
	{
	  err.h = h[j];
	  double steps = (xMax - xMin)/h[j];
	  err.xVals.clear();
	  err.derExact.clear();
	  err.derApprox.clear();
	  err.errorNorms.clear();
	  err.errorsX.clear();
	  for(int k = 0; k <= steps; ++k)
	    {
	      double x = xMin + h[j] * k;
	      double exact = 0;
	      double approx = 0;
	      double errorNorm = 0;
	      double error = 0;
	      if(i == 0) // CD
		{
		  exact = fn_prime(x);
		  approx = CentralDifference(x, h[j], fn);
		  errorNorm = RelError(exact, approx);
		  error = Error(exact, approx);
		}
	      if(i == 1) // FWD
		{
		  exact = fn_prime(x);
		  approx = ForwardDifference(x, h[j], fn);
		  errorNorm = RelError(exact, approx);
		  error = Error(exact, approx);
		}
	      err.xVals.push_back(x);
	      err.derExact.push_back(exact);
	      err.derApprox.push_back(approx);
	      err.errorNorms.push_back(errorNorm);
	      err.errorsX.push_back(error);
	    }
	  WriteErrors(err);
	  cout << "Error Norm fn'(x) " << err.derLabel << " [h = " << err.h << "]: " << ErrorNorm(err.derExact, err.derApprox) << endl;
	  vector<double>::iterator result = max_element(err.errorsX.begin(), err.errorsX.end());
	  int index = distance(err.errorsX.begin(), result);
	  cout << '\t' << "Max Error: " << err.errorsX[index] << endl;
	}
    }
  // fn''
  for(int i = 0; i < h.size(); ++i)
    {
      err.derLabel = labels[2];
      err.derType = "fnPP";
      double steps = (xMax - xMin)/h[i];
      err.h = h[i];
      err.xVals.clear();
      err.derExact.clear();
      err.derApprox.clear();
      err.errorNorms.clear();
      err.errorsX.clear();
      for(int j = 0; j <= steps; ++j)
	{
	  double x = xMin + h[i] * j;
	  double exact = fn_pprime(x);
	  double approx = CentralDifference_2(x, h[i], fn);
	  double errorNorm = RelError(exact, approx);
	  double error = Error(exact, approx);
	  err.xVals.push_back(x);
	  err.derExact.push_back(exact);
	  err.derApprox.push_back(approx);
	  err.errorNorms.push_back(errorNorm);
	  err.errorsX.push_back(error);
	}
      WriteErrors(err);
      cout << "Error Norm fn''(x) " << err.derLabel << " [h = " << err.h << "]: " << ErrorNorm(err.derExact, err.derApprox) << endl;
      vector<double>::iterator result = max_element(err.errorsX.begin(), err.errorsX.end());
      int index = distance(err.errorsX.begin(), result);
      cout << '\t' << "Max Error: " << err.errorsX[index] << endl;
    }
  return 0;
}

double fn(double x)
{
  return sin(2*M_PI*x)+5;
}

double fn_prime(double x)
{
  return 2*M_PI*cos(2*M_PI*x);
}

double fn_pprime(double x)
{
  return -4*M_PI*M_PI*sin(2*M_PI*x);
}

double RelError(double exact, double approx)
{
  /*
  double error = abs((approx - exact)/approx);
  if(error > .9)
    error = 0.0;
  return error;
  */
  return abs((approx - exact)/approx);
}

double Error(double exact, double approx)
{
  return abs(approx-exact);
}

double ErrorNorm(vector<double> exact, vector<double> approx)
{
  double sum = 0.0;
  for(int i = 0; i < exact.size(); ++i)
    {
      double error = Error(exact[i], approx[i]);
      sum += pow(error, 2);
    }
  return pow((1.0/exact.size())*sum, 1.0/2.0);
}

template <typename T>
void WriteFile(vector<T> x, vector<T> y, T h, string filename, bool exact)
{
  ostringstream hSS;
  string hStr;
  string outFile = "output_data_";
  if(exact == false)
    {
      hSS << h;
      hStr = hSS.str();
      outFile += filename + "_h_" + hStr + ".dat";
    }
  else
    outFile += filename + ".dat";
    
  ofstream output;
  output.open(outFile.c_str());
  for(int i = 0; i < x.size(); ++i)
    {
      output << x[i] << '\t' << y[i] << endl;
    }
  output.close();
}

template <typename T>
void WriteErrors(errors<T> err)
{
  ostringstream hSS;
  string hStr;
  string outFile = "output_error_";
  hSS << err.h;
  hStr = hSS.str();
  outFile += err.derType + "_" + err.derLabel + "_h_" + hStr + ".dat";
  ofstream output;
  output.open(outFile.c_str());
  for(int i = 0; i < err.xVals.size(); ++i)
    {
      T x = err.xVals[i];
      T exact = err.derExact[i];
      T approx = err.derApprox[i];
      T errorNorm = err.errorNorms[i];
      T error = err.errorsX[i];
      output << x << '\t' << exact << '\t' << approx << '\t' << errorNorm << '\t' << error << endl;
    }
  output.close();
}
