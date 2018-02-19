#include <iostream>
#include <string>
#include "cubic_spline.h"
#include "misc.h"
#include "io.h"

int main(int argc, char** argv){
  Real x[4]={1,2,3,4};
  Real y[4]={1,8,27,64};
  CubicSpline<Real> cs(x,y,4);
  Real x_test=1.2;
  std::cout<<"Cubic spline value: "<<cs(x_test)<<std::endl;
  std::cout<<"Accurate value: "<<x_test*x_test*x_test<<std::endl;
  x_test=5;
  std::cout<<"Cubic spline value: "<<cs(x_test)<<std::endl;
  std::cout<<"Accurate value: "<<x_test*x_test*x_test<<std::endl;

  std::string fname(argv[1]);
  std::vector<Real> re_sts_mat;
  int n_row,n_col;
  char sep=' ';
  int skiprows=14;
  ReadMatrix(fname,re_sts_mat,n_row,n_col,sep,skiprows);
  Real* yplus_arr=new Real[n_row];
  Real* uvplus_arr=new Real[n_row];
  for(int i=0;i<n_row;i++){
    yplus_arr[i]=re_sts_mat[i*n_col+1];
    uvplus_arr[i]=re_sts_mat[i*n_col+6];
  }
  CubicSpline<Real> cs2(yplus_arr,uvplus_arr,n_row);
  Real yplus_test=79.0;
  std::cout<<"Cubic spline value: "<<cs2(yplus_test)<<std::endl;
  delete[] yplus_arr;
  delete[] uvplus_arr;
  return 0;
}
