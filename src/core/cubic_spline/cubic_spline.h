#ifndef CUBIC_SPLINE
#define CUBIC_SPLINE

template<typename T>
class CubicSpline {
public:
  CubicSpline(){}

  CubicSpline(T* in_x_arr,T* in_y_arr,int in_n,int in_bc_type=0,
    int in_extrapolate_type=0);

  ~CubicSpline();

  void Init(T* in_x_arr,T* in_y_arr,int in_n,int in_bc_type,
    int in_extrapolate_type);

  void SetCoefMat(int in_bc_type);

  void SolveCoef();

  T operator()(T x);

private:
  int n;
  T* x_arr=nullptr;
  T* y_arr=nullptr;
  T* A_mat=nullptr;
  T* b_arr=nullptr;
  T* y_deriv2_arr=nullptr;
  T* dx_arr=nullptr;
  T* dydx_arr=nullptr;
  int bc_type;
  int extrapolate_type;
};
#endif