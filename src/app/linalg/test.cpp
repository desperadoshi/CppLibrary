#include <iostream>
#include <string>
#include "misc.h"
#include "linalg.h"

int main(int argc, char** argv){
  /*
   *
   */
  Real a=1.0;
  std::cout<<sgn(a)<<std::endl;
  int n=4;
  Real* x_arr=new Real[n];
  linspace(x_arr,0.0,1.0,n);
  for(int i=0;i<n;i++){
    std::cout<<x_arr[i]<<",";
  }
  std::cout<<std::endl;
  logspace(x_arr,0.0,1.0,n);
  for(int i=0;i<n;i++){
    std::cout<<x_arr[i]<<",";
  }
  std::cout<<std::endl;
  // error_here();
  /*
   *
   */
  delete[] x_arr;
  return 0;
}
