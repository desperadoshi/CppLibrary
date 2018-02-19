#include <iostream>
#include <string>
#include "misc.h"
#include "linalg.h"
#include "io.h"

int main(int argc, char** argv){
  int n=4;
  Real* x_arr=new Real[n];
  logspace(x_arr,0.0,1.0,n);
  /*
   *
   */
  std::string fname="test.dat";
  char sep=' ';
  std::string header="# This is a line of the header";
  WriteMatrix(fname,x_arr,2,2,9,sep,header);
  std::vector<Real> x_new_arr;
  int n_row,n_col;
  int skiprows=1;
  ReadMatrix(fname,x_new_arr,n_row,n_col,sep,skiprows);
  for(auto x : x_new_arr){
    std::cout<<x<<sep;
  }
  std::cout<<std::endl;

  WriteMatrix(fname,x_arr,2,2);
  x_new_arr.clear();
  ReadMatrix(fname,x_new_arr,n_row,n_col);
  for(auto x : x_new_arr){
    std::cout<<x<<sep;
  }
  std::cout<<std::endl;
  /*
   *
   */
  delete[] x_arr;
  return 0;
}
