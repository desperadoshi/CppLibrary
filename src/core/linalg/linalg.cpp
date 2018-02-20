#ifndef LINALG_H
#define LINALG_H

#include "misc.h"
#include <cmath>
#include <cassert>
#include <limits>

template<typename T>
void linspace(T* lin_arr,T min,T max,int n_pts){
    assert(n_pts>1);
    T step=(max-min)/(n_pts-1);
    for(int i=0;i<n_pts;i++){
        lin_arr[i]=min+i*step;
    }
} // linspace

template<typename T>
void logspace(T* log_arr,T min,T max,int n_pts){
    assert(n_pts>1);
    T step=(max-min)/(n_pts-1);
    for(int i=0;i<n_pts;i++){
        log_arr[i]=pow(10,min+i*step);
    }
} // logspace

/** \brief LU decomposition. Note: This routine is not reliable.
 * @param a    The matrix to be decomposed. An array of shape (n,n).
 * @param n    Indicate the size
 * @param indx An array of shape (n). The output vector that records the row
 *             permutation effected by the parital pivoting.
 * @param d    +1 or -1, indicating whether the number of row interchanges was
 *             even or odd.
 *
 * From Numerical Recipe.
 * Cubic spline tests show that this routine along with LUbksb is not reliable.
 * Jingchang Shi 2018-02-20, 10:40:05
 */
template<typename T>
void LUdcmp(T** a,int n,int* indx,T* d){
  int i,imax,j,k;
  T big=0.0,dum=0.0,sum=0.0;
  T *vv = new T[n];

  *d = 1.;
  for(i=0;i<n;i++) {
    big=0.;
    for(j=0;j<n;j++)
      //if( (temp=fabs(a[i][j])) > big) big=temp;
      if(fabs(a[i][j]) > big)
        big=fabs(a[i][j]);
    if(big == 0.) {
      error_here("Singular Matrix!");
    }
    vv[i]=1./big;
  }

  for(j=0;j<n;j++) {
    for(i=0;i<j;i++) {
      sum=a[i][j];
      for(k=0;k<i;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.;
    for(i=j;i<n;i++) {
      sum=a[i][j];
      for(k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if((dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }

    if(j != imax) {
      for(k=0;k<n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if(a[j][j] == 0.0)
      a[j][j]=std::numeric_limits<T>::min();
    if(j != n-1) {
      dum=1./a[j][j];
      for(i=j+1;i<n;i++)
        a[i][j] *= dum;
    }
  }

  delete [] vv;
} // LUdcmp

/** \brief Solves Ax=b. Note: This routine is not reliable.
 * @param a    LU decomposition of a matrix from ludcmp
 * @param indx A vector, created by LUdcmp, containing the row permutations
 *             effected by the partial pivoting.
 * @param b    rhs. Also the output. b is replaced by the solution vector.
 *
 * Cubic spline tests show that this routine along with LUdcmp is not reliable.
 * The solution from this routine is different from what is obtained by Eigen3.
 * Jingchang Shi 2018-02-20, 10:40:05
 */
template<typename T>
void LUbksb(T** a,int n,int* indx,T* b){
  int i,ii=-1,ip,j;
  T sum;

  // std::ios::fmtflags old_settings = std::cout.flags();
  // int old_precision = std::cout.precision();
  // std::cout<<setprecision(std::numeric_limits<T>::digits10+2)<<scientific;

  for(i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if(ii>-1)
      for(j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if(sum) ii=i;
    b[i]=sum;
  }

  for(i=n-1;i>=0;i--) {
    sum=b[i];
    for(j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }

  // std::cout.flags(old_settings);
  // std::cout.precision(old_precision);
} // LUbksb

template void linspace(Real* lin_arr,Real min,Real max,int n_pts);
template void logspace(Real* log_arr,Real min,Real max,int n_pts);
template void LUdcmp(Real** a,int n,int* indx,Real *d);
template void LUbksb(Real** a,int n,int* indx,Real *b);

#endif
