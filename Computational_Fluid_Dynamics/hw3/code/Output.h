#ifndef OUTPUT_H
#define OUTPUT_H
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "Dimensions.h"
#include "Matrix.h"

void WriteFile(Dimensions<double> &pos, std::vector<double> &omega, std::vector<double> &psi, std::vector<double> u, std::vector<double> v, std::string filename)
{
  std::stringstream ss;
  //string outputfile = filename;
  std::ofstream outfile;
  outfile.open(filename.c_str());
  outfile << std::setprecision(6) << std::fixed;
  
  for(int i = 0; i < pos.y.size(); ++i)
    {
      for(int j = 0; j < pos.x.size(); ++j)
	outfile << pos.x[j] << ',' << pos.y[i] << ',' << omega[j + i * pos.x.size()] << ',' << psi[j + i * pos.x.size()] << ',' << u[j + i * pos.x.size()] << ',' << v[j + i * pos.x.size()] << std::endl;
    }
  outfile.close();
  return;
}

void WriteVorticity(Dimensions<double> &pos, std::vector<double> &omega, std::string filename)
{
  std::ofstream outfile;
  outfile.open(filename.c_str());
  outfile << std::setprecision(6) << std::fixed;

  for(int i = 0; i < pos.y.size(); ++i)
    {
      for(int j = 0; j < pos.x.size(); ++j)
	outfile << pos.x[j] << ',' << pos.y[i] << ',' << omega[j + i * pos.x.size()] << std::endl;
    }
  outfile.close();
}

template <typename T>
void WriteMatrix(Matrix<T> &A, const int &Nx)
{
  // Print the matrix in latex format
  std::string filename = "matrix.tex";
  std::string cmd = "pdflatex ";
  cmd.append(filename);
  std::ofstream output;
  
  output.open(filename.c_str());

  output << "\\documentclass{article}" << std::endl;
  output << "\\usepackage{amsmath}" << std::endl;
  output << "\\usepackage[a4paper, landscape, margin=1in]{geometry}" << std::endl;
  output << "\\begin{document}" << std::endl << std::endl;
  output << "\\thispagestyle{empty}" << std::endl;
  output << "\\[" << std::endl;
  output << "\\mathbf{A} = \\left[ \\begin{array}{";
  for(int i = 0; i < Nx*Nx; ++i)
    output << "c ";
  output << "}" << std::endl;

  for(int i = 0; i < Nx*Nx; ++i)
    {
      std::string line = "";
      int counter = 0;

      // Lower Lower Diagonal
      if(i >= Nx)
	{
	  if(i > Nx)
	    {
	      for(int j = 0; j < i - Nx; ++j)
		{
		  line.append(" 0 & ");
		  ++counter;
		}
	    }
	  if((int)A.llower[i] >= 0)
	    line.append(" ");
	  line.append(std::to_string((int)A.llower[i]));
	  line.append(" & ");
	  ++counter;
	  for(int j = 0; j < Nx - 2; ++j)
	    {
	      line.append(" 0 & ");
	      ++counter;
	    }
	}
      else
	{
	  for(int j = 0; j < i - 1; ++j)
	    {
	      line.append(" 0 & ");
	      ++counter;
	    }
	}

      // Lower Diagonal
      if(i > 0)
	{
	  if((int)A.lower[i] >= 0)
	    line.append(" ");
	  line.append(std::to_string((int)A.lower[i]));
	  line.append(" & ");
	  ++counter;
	}

      // Diagonal
      line.append(" ");
      line.append(std::to_string((int)A.diag[i]));
      if(i < Nx*Nx - 1)
	line.append(" & ");
      ++counter;

      // Upper Diagonal
      if(i < Nx*Nx - 1)
	{
	  if((int)A.upper[i] >= 0)
	    line.append(" ");
	  line.append(std::to_string((int)A.upper[i]));
	  if(i < Nx*Nx - 2)
	    line.append(" & ");
	  ++counter;
	}

      // Upper Upper Diagonal
      if(i < Nx*Nx - Nx)
	{
	  for(int j = 0; j < Nx - 2; ++j)
	    {
	      line.append(" 0 & ");
	      ++counter;
	    }
	  if((int)A.uupper[i] >= 0)
	    line.append(" ");
	  line.append(std::to_string((int)A.uupper[i]));
	  if(i < Nx*Nx - Nx - 1)
	    line.append(" & ");
	  ++counter;
	}

      for(int j = counter; j < Nx*Nx - 1; ++j)
	{
	  line.append(" 0 & ");
	  ++counter;
	}

      if(i < Nx*Nx - 2 && i != Nx*Nx - Nx - 1)
	line.append(" 0");

      if(i < Nx*Nx - 1)
	line.append("\\\\");
     
      output << line << std::endl;
    }
  output << "\\end{array}\\right]" << std::endl;
  output << "\\]" << std::endl << std::endl;
  output << "\\end{document}" << std::endl;
  output.close();

  system(cmd.c_str());

  return;
}


#endif
