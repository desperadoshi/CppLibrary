#include <iostream>
#include "cubic_spline.h"

typedef double Real;

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
  return 0;
}