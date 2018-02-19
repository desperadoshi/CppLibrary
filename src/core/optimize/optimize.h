#ifndef OPTIMIZE_H
#define OPTIMIZE_H
/*
 * Not Tested.
 */
/*
 * Newton method to find the root of nonlinear equation.
 * The derivative is approximated by a central difference.
 */
template<typename T>
T Root1DNewton(T func(T x, const T* pmt),
    const T x0, const T tol, const T* pmt){
    T root=x0;
    const int MAXIT=20;
    const T step=tol;
    for(int i=0;i<MAXIT;i++){
        T dfdx=(func(root+step,pmt)-func(root-step,pmt))/(2.0*step);
        assert(dfdx!=0);
        T func_val=func(root,pmt);
        T diff=func_val/dfdx;
        // cout<<diff<<endl;
        root-=diff;
        // if(fabs(diff)<tol){
        if(fabs(func_val)<tol){
            cout<<i<<endl;
            return root;
        }
    }
    cout<<"In the function Root1DNewton, MAX iteration is reached."<<endl;
    exit(1);
} // Root1DNewton

template<typename T>
void Root2DNewton(void func(T* XVec,T* YVec,const T* pmt),
    T* root,const T tol,const T* pmt){
    const int MAXIT=100;
    const T step=tol;
    T XVec1[2],YVec1[2],XVec2[2],YVec2[2];
    T diff[2];
    T a11,a12,a21,a22;
    for(int i=0;i<MAXIT;i++){
        XVec1[0]=root[0]-step;
        XVec1[1]=root[1];
        func(XVec1,YVec1,pmt);
        XVec2[0]=root[0]+step;
        XVec2[1]=root[1];
        func(XVec2,YVec2,pmt);
        a11=(YVec2[0]-YVec1[0])/(2.0*step);
        a21=(YVec2[1]-YVec1[1])/(2.0*step);
        XVec1[0]=root[0];
        XVec1[1]=root[1]-step;
        func(XVec1,YVec1,pmt);
        XVec2[0]=root[0];
        XVec2[1]=root[1]+step;
        func(XVec2,YVec2,pmt);
        a12=(YVec2[0]-YVec1[0])/(2.0*step);
        a22=(YVec2[1]-YVec1[1])/(2.0*step);
        func(root,YVec1,pmt);
        diff[0]=(a12*YVec1[1]-a22*YVec1[0])/(a11*a22-a12*a21);
        diff[1]=(a21*YVec1[0]-a11*YVec1[1])/(a11*a22-a12*a21);
        root[0]+=diff[0];
        root[1]+=diff[1];
        cout<<i<<endl;
        cout<<diff[0]<<", "<<diff[1]<<endl;
        if(fabs(diff[0])<tol and fabs(diff[1])<tol){
            return;
        }
    }
    cout<<"In the function Root1DNewton, MAX iteration is reached."<<endl;
    exit(1);
} // Root2DNewton

Real ZRIDDR(Real (*func)(const Real, Real*), const Real eq_lhs,
            const Real x1, const Real x2, const Real xacc, Real* pmt) {
    const int MAXIT = 60;
    Real fl = func(x1, pmt)-eq_lhs;
    Real fh = func(x2, pmt)-eq_lhs;
    if((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
        Real xl = x1;
        Real xh = x2;
        Real ans = -9.99E99;
        for(int j=0; j<MAXIT; j++) {
            Real xm = 0.5 * (xl + xh);
            Real fm = func(xm, pmt)-eq_lhs;
            Real s = sqrt(fm*fm-fl*fh);
            if(s == 0.0) {
                // cout<<j+1<<endl;
                return ans;
            }
            Real xnew = xm + (xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
            // no matching function for call to ‘abs(hpMusic::Dual)’
            if(fabs(xnew-ans) <= xacc) {
                // cout<<j+1<<endl;
                return ans;
            }
            ans=xnew;
            Real fnew = func(ans, pmt)-eq_lhs;
            if(fnew==0.0) {
                // cout<<j+1<<endl;
                return ans;
            }
            if(copysign(fm, fnew) != fm) {
                xl=xm;
                fl=fm;
                xh=ans;
                fh=fnew;
            } else if(copysign(fl, fnew) != fl) {
                xh=ans;
                fh=fnew;
            } else if(copysign(fh, fnew) != fh) {
                xl=ans;
                fl=fnew;
            } else throw("never get here.");
            // no matching function for call to ‘abs(hpMusic::Dual)’
            if(fabs(xh-xl) <= xacc) {
                // cout<<j+1<<endl;
                return ans;
            }
        }
        throw("zriddr exceed maximum iterations.");
    } else {
        if(fl == 0.0) return x1;
        if(fh == 0.0) return x2;
        cout<<"fl("<<x1<<")="<<fl<<endl;
        cout<<"fh("<<x2<<")="<<fh<<endl;
        throw("root must be bracketed in zriddr.");
    }
} // ZRIDDR

#endif OPTIMIZE_H