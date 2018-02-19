#ifndef LINALG_H
#define LINALG_H

#include "misc.h"
#include <cmath>
#include <cassert>

template<typename T>
T linspace(T* lin_arr,T min,T max,int n_pts){
    assert(n_pts>1);
    T step=(max-min)/(n_pts-1);
    for(int i=0;i<n_pts;i++){
        lin_arr[i]=min+i*step;
    }
} // linspace

template<typename T>
T logspace(T* log_arr,T min,T max,int n_pts){
    assert(n_pts>1);
    T step=(max-min)/(n_pts-1);
    for(int i=0;i<n_pts;i++){
        log_arr[i]=pow(10,min+i*step);
    }
} // logspace

template Real linspace(Real* lin_arr,Real min,Real max,int n_pts);
template Real logspace(Real* log_arr,Real min,Real max,int n_pts);

#endif