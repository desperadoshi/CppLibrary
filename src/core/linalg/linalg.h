#ifndef LINALG_H
#define LINALG_H

template <typename T>
inline int sgn(T val) {
    return (val>T(0))-(val<T(0));
} // sgn

template<typename T>
void linspace(T* lin_arr,T min,T max,int n_pts);

template<typename T>
void logspace(T* log_arr,T min,T max,int n_pts);

template<typename T>
void LUdcmp(T** a,int n,int* indx,T* d);

template<typename T>
void LUbksb(T** a,int n,int* indx,T* b);

#endif
