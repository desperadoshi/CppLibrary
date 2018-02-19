#include <iostream>
#include <string>
#include "misc.h"
#include "linalg.h"
#include "io.h"

int main(int argc, char** argv){
  Real a=1.0;
  std::cout<<sgn(a)<<std::endl;
  int n=4;
  Real* x_arr=new Real[n];
  logspace(x_arr,0.0,1.0,n);
  /*
   *
   */
  std::string fname="test.dat";
  char sep=' ';
  WriteMatrix(fname,x_arr,2,2,9,sep);
  std::vector<Real> x_new_arr;
  ReadMatrix(fname,x_new_arr,2,2,sep);
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