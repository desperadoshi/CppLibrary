#ifndef LINALG_H
#define LINALG_H

template <typename T>
inline int sgn(T val) {
    return (val>T(0))-(val<T(0));
} // sgn

template<typename T>
T linspace(T* lin_arr, T min, T max, int n_pts);

template<typename T>
T logspace(T* log_arr, T min, T max, int n_pts);

#endif