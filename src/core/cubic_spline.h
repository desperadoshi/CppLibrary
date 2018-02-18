#ifndef CUBIC_SPLINE
#define CUBIC_SPLINE

#include <algorithm>
#include <iostream>
#include <cassert>

#include <Eigen/Dense>

template<typename T>
class CubicSpline {
public:
  CubicSpline(){}

  CubicSpline(T* in_x_arr,T* in_y_arr,int in_n,int in_bc_type=0,
    int in_extrapolate_type=0){
    Init(in_x_arr,in_y_arr,in_n,in_bc_type,in_extrapolate_type);
    SetCoefMat(in_bc_type);
    SolveCoef();
  }

  ~CubicSpline(){
    if(x_arr){
      delete[] x_arr;
      x_arr=NULL;
    }
    if(y_arr){
      delete[] y_arr;
      y_arr=NULL;
    }
    if(A_mat){
      delete[] A_mat;
    }
    if(b_arr){
      delete[] b_arr;
    }
    if(y_deriv2_arr){
      delete[] y_deriv2_arr;
    }
    if(dx_arr){
      delete[] dx_arr;
    }
    if(dydx_arr){
      delete[] dydx_arr;
    }
  }

  void Init(T* in_x_arr,T* in_y_arr,int in_n,int in_bc_type,int in_extrapolate_type){
    n=in_n;
    x_arr=new T[n];
    y_arr=new T[n];
    int* idx=new int[n];
    for(int i=0;i<n;i++){
      idx[i]=i;
    }
    std::sort(idx,idx+n,[&in_x_arr](int i1,int i2){return in_x_arr[i1]<in_x_arr[i2];});
    for(int i=0;i<n;i++){
      x_arr[i]=in_x_arr[idx[i]];
      y_arr[i]=in_y_arr[idx[i]];
    }
    delete[] idx;
    A_mat=new T[n*n];
    b_arr=new T[n];
    y_deriv2_arr=new T[n];
    dx_arr=new T[n-1];
    dydx_arr=new T[n-1];
    bc_type=in_bc_type;
    extrapolate_type=in_extrapolate_type;
  }

  void SetCoefMat(int in_bc_type){
    for(int i=0;i<n-1;i++){
      dx_arr[i]=x_arr[i+1]-x_arr[i];
      dydx_arr[i]=(y_arr[i+1]-y_arr[i])/dx_arr[i];
    }
    for(int i=1;i<n-1;i++){
      A_mat[i*n+i-1]=dx_arr[i-1]/6.0;
      A_mat[i*n+i]=(dx_arr[i]+dx_arr[i-1])/3.0;
      A_mat[i*n+i+1]=dx_arr[i]/6.0;
      b_arr[i]=dydx_arr[i]-dydx_arr[i-1];
    }
    /*
     * Not-A-knot
     * Set the 3rd deriv to be continuous.
     */
    if(in_bc_type==0){
      A_mat[0]=dx_arr[1];
      A_mat[1]=-(dx_arr[1]+dx_arr[0]);
      A_mat[2]=dx_arr[0];
      b_arr[0]=0.0;
      A_mat[n*n-3]=dx_arr[n-2];
      A_mat[n*n-2]=-(dx_arr[n-2]+dx_arr[n-3]);
      A_mat[n*n-1]=dx_arr[n-3];
      b_arr[n-1]=0.0;
    /*
     * Set the value of the 1st deriv
     */
    } else if(in_bc_type==1){
      std::cerr<<"BC of the 1st deriv is not implemented!"<<std::endl;
    /*
     * Set the value of the 2nd deriv
     */
    } else if(in_bc_type==2){
      std::cerr<<"BC of the 2nd deriv is not implemented!"<<std::endl;
    }
  }

  void SolveCoef(){
    Eigen::MatrixXd A(n,n);
    Eigen::VectorXd b(n);
    for(int i=0;i<n;i++){
      b(i)=b_arr[i];
      for(int j=0;j<n;j++){
        A(i,j)=A_mat[i*n+j];
      }
    }
    Eigen::VectorXd X(n);
    X=A.colPivHouseholderQr().solve(b);
    for(int i=0;i<n;i++){
      y_deriv2_arr[i]=X(i);
    }
  }

  T operator()(T x){
    T y;
    if(x<x_arr[0] or x>x_arr[n-1]){
      if(extrapolate_type==0){
        y=(x<x_arr[0])*y_arr[0]+(x>x_arr[n-1])*y_arr[n-1];
        return y;
      }
    }
    int j=std::lower_bound(x_arr,x_arr+n,x)-x_arr;
    j-=1;
    T A=(x_arr[j+1]-x)/dx_arr[j],B=1-A;
    T C=(A*A*A-A)*dx_arr[j]*dx_arr[j]/6.0;
    T D=(B*B*B-B)*dx_arr[j]*dx_arr[j]/6.0;
    y=A*y_arr[j]+B*y_arr[j+1]+C*y_deriv2_arr[j]+D*y_deriv2_arr[j+1];
    return y;
  }

private:
  int n;
  T* x_arr=NULL;
  T* y_arr=NULL;
  T* A_mat=NULL;
  T* b_arr=NULL;
  T* y_deriv2_arr=NULL;
  T* dx_arr=NULL;
  T* dydx_arr=NULL;
  int bc_type;
  int extrapolate_type;
};
#endif