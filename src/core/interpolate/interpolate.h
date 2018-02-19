#ifndef INTERPOLATE_H
#define INTERPOLATE_H
/*
 * Not Tested.
 */
/**
 * Return the coefficients of 1D polynomial
 * Vieta's formula relates the roots with the coefficients of the polynomial.
 * According to stackoverflow, employing FFT could speed up to O(nlogn).
 */
template<typename T>
void GetPoly1DCoefByRoots(const T* root_arr, const int N,
    Real* coef_arr){
    coef_arr[0]=-root_arr[0];
    coef_arr[1]=1;
    for(int k=2;k<=N;k++){
        coef_arr[k]=1;
        for(int i=k-2;i>=0;i--){
            coef_arr[i+1]=coef_arr[i]-root_arr[k-1]*coef_arr[i+1];
        }
        coef_arr[0]*=-root_arr[k-1];
    }
} // GetPoly1DCoefByRoots

/**
 * Return the coefficients of the Lagrange 1D polynomial
 * The order is N-1, where N is the number of points
 */
template<typename T>
void Lagrange1DBasisPolyCoef(const T* x_arr, const int nx, const int j,
    Real* coef_arr){
    T root_arr[nx-1];
    T comm=1;
    int count=0;
    for(int i=0;i<nx;i++){
        if(i==j) continue;
        root_arr[count++]=x_arr[i];
        comm*=x_arr[j]-x_arr[i];
    }
    assert(count==nx-1);
    GetPoly1DCoefByRoots(root_arr,nx-1,coef_arr);
    for(int i=0;i<nx;i++){
        coef_arr[i]/=comm;
    }
} // Lagrange1DBasisPolyCoef

template<typename T>
void Lagrange1DBasisDerivPolyCoef(const T* x_arr, const int nx, const int j,
    Real* deriv_coef_arr){
    Real coef_arr[nx];
    Lagrange1DBasisPolyCoef(x_arr,nx,j,coef_arr);
    for(int i=1;i<nx;i++){
        deriv_coef_arr[i-1]=coef_arr[i]*i;
    }
} // Lagrange1DBasisDerivPolyCoef

template<typename T>
void Lagrange1DDerivPolyCoef(const T* x_arr, const T* y_arr, const int nx,
    Real* deriv_coef_arr){
    for(int i=0;i<nx-1;i++){
        deriv_coef_arr[i]=0.0;
    }
    Real basis_deriv_coef_arr[nx-1];
    for(int j=0;j<nx;j++){
        Lagrange1DBasisDerivPolyCoef(x_arr,nx,j,basis_deriv_coef_arr);
        for(int i=0;i<nx-1;i++){
            deriv_coef_arr[i]+=basis_deriv_coef_arr[i]*y_arr[j];
        }
    }
} // Lagrange1DDerivPolyCoef

template<typename T>
T Interp1DByBasis(const T x, const T* coef, const int nx, const T* y){
    /*
     * coef has nx*nx elements.
     * y has nx elements.
     */
    T sum=0.0;
    for(int i=0;i<nx;i++){
        for(int j=0;j<nx;j++){
            sum+=y[i]*(coef[i*nx+j]*pow(x,j));
        }
    }
    return sum;
} // Interp1DByBasis

/**
 * Make sure arr_x is sorted
 */
template<typename T>
T Interp1D(T* arr_x, T* arr_y, T x0, int npts){
    T* low=lower_bound(&arr_x[0],&arr_x[npts-1],x0);
    T* up=upper_bound(&arr_x[0],&arr_x[npts-1],x0);
    T* pos_l=low==up?low-1:low;
    T* pos_h=low==up?low:up;
    // linear interpolation
    int idx_l=pos_l-&arr_x[0];
    int idx_h=pos_h-&arr_x[0];
    T dx=arr_x[idx_h]-arr_x[idx_l];
    assert(dx!=0.0);
    T y0=arr_y[idx_l]+(arr_y[idx_h]-arr_y[idx_l])/dx*(x0-arr_x[idx_l]);
    return y0;
}

#endif INTERPOLATE_H