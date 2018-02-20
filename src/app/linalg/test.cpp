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
  delete[] x_arr;
  /*
   *
   */
  Real** A_mat=new Real*[2];
  A_mat[0]=new Real[2];
  A_mat[1]=new Real[2];
  A_mat[0][0]=1.0;
  A_mat[0][1]=3.0;
  A_mat[1][0]=3.0;
  A_mat[1][1]=2.0;
  int indx[2];
  Real d;
  LUdcmp(A_mat,2,indx,&d);
  std::cout<<A_mat[0][0]<<","<<A_mat[0][1]<<std::endl;
  std::cout<<A_mat[1][0]<<","<<A_mat[1][1]<<std::endl;
  std::cout<<indx[0]<<","<<indx[1]<<std::endl;
  std::cout<<d<<std::endl;
  Real b_arr[2]={1.0,0.9};
  LUbksb(A_mat,2,indx,b_arr);
  std::cout<<b_arr[0]<<","<<b_arr[1]<<std::endl;
  delete[] A_mat[0];
  delete[] A_mat[1];
  delete[] A_mat;
  return 0;
}
