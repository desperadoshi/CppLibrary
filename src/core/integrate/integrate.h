#ifndef INTEGRATE_H
#define INTEGRATE_H
/*
 * Not Tested.
 */

void RK45(Real (*func)(const Real,const Real,const Real*),
          Real& x,Real& y,const Real h,const Real* pmt) {
    Real k1=h*func(x,y,pmt);
    Real k2=h*func(x+0.5*h,y+0.5*k1,pmt);
    Real k3=h*func(x+0.5*h,y+0.5*k2,pmt);
    Real k4=h*func(x+h,y+k3,pmt);
    y+=k1/6.0+k2/3.0+k3/3.0+k4/6.0;
    x+=h;
} // RK45

void RKAdaptive(Real (*func)(const Real,const Real,const Real*),
    const Real x0,const Real xend,Real& y,Real h,const Real acc,
    const Real* pmt) {
    // Real y0=y;
    Real x1,x2;
    Real y1,y2;
    Real x=x0;
    Real err;
    int count;
    const int MAXIT=20;
    Real h_factor=1.0;
    Real scale=1.0;
    while(fabs(x-xend)/xend>1E-9) {
    // while(x<xend) {
        // Take full step
        x1=x;
        y1=y;
        RK45(func,x1,y1,h,pmt);
        // Take two half-steps
        x2=x;
        y2=y;
        RK45(func,x2,y2,h/2,pmt);
        RK45(func,x2,y2,h/2,pmt);
        /*
         * Use absolute truncation error
         */
        err=fabs(y2-y1)/15;
        /*
         * Use relative error
         */
        // err=fabs((y2-y1)/y2);
        count=0;
        // cout<<setprecision(15)<<xend<<endl;
        // cout<<setprecision(15)<<x<<endl;
        while(err>acc) {
            // cout<<setprecision(15)<<err<<endl;
            // cout<<setprecision(15)<<h<<endl;
            count++;
            if(count>MAXIT) {
                cout<<"In RKAdaptive, MAX iteration is reached."<<endl;
                cout<<"Current error is "<<err<<", the goal is "<<acc<<endl;
                // break;
                exit(1);
            }
            // Calculate new step-length
            // h_factor=pow(0.2,tanh(err/acc-1.0));
            // h*=pow(fabs(acc/err),h_factor)*scale;
            h*=scale*min(max(acc/err,0.3),2.0);
            h=(x+h)>xend?(xend-x):h;
            // Take full step
            x1=x;
            y1=y;
            RK45(func,x1,y1,h,pmt);
            // Take two half-steps
            x2=x;
            y2=y;
            RK45(func,x2,y2,h/2,pmt);
            RK45(func,x2,y2,h/2,pmt);
            /*
             * Use absolute truncation error
             */
            err=fabs(y2-y1)/15;
            /*
             * Use relative error
             */
            // err=fabs((y2-y1)/y2);
        }
        // cout<<count<<endl;
        y=y2;
        x=x2;
        // Calculate new step-length
        // h_factor=pow(0.2,tanh(err/acc-1.0));
        // h*=pow(acc/err,h_factor)*scale;
        h*=scale*min(max(acc/err,0.3),2.0);
        h=(x+h)>xend?(xend-x):h;
    }
} // RKAdaptive

/*
 * Cash-Karp Runge-Kutta, 5th order.
 *
 */
void RKCK(Real (*func)(const Real,const Real,const Real*),
          Real x,Real y,const Real h,const Real* pmt,Real& ytemp,Real& yerr){
    const Real a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,
               b21=0.2,
               b31=3.0/40.0,b32=9.0/40.0,
               b41=0.3,b42=-0.9,b43=1.2,
               b51=-11.0/54.0,b52=2.5,b53=-70.0/27.0,b54=35.0/27.0,
               b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,b64=44275.0/110592.0,b65=253.0/4096.0,
               c1=37.0/378.0,c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0;
    const Real dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,dc4=c4-13525.0/55296.0,
               dc5=-277.00/14336.0,dc6=c6-0.25;
    Real ak1,ak2,ak3,ak4,ak5,ak6;
    ak1=func(x,y,pmt);
    ytemp=y+h*b21*ak1;
    ak2=func(x+a2*h,ytemp,pmt);
    ytemp=y+h*(b31*ak1+b32*ak2);
    ak3=func(x+a3*h,ytemp,pmt);
    ytemp=y+h*(b41*ak1+b42*ak2+b43*ak3);
    ak4=func(x+a4*h,ytemp,pmt);
    ytemp=y+h*(b51*ak1+b52*ak2+b53*ak3+b54*ak4);
    ak5=func(x+a5*h,ytemp,pmt);
    ytemp=y+h*(b61*ak1+b62*ak2+b63*ak3+b64*ak4+b65*ak5);
    ak6=func(x+a6*h,ytemp,pmt);
    ytemp=y+h*(c1*ak1+c3*ak3+c4*ak4+c6*ak6);
    yerr=h*(dc1*ak1+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6);
} // RKCK

void RKQS(Real (*func)(const Real,const Real,const Real*),
    const Real x0,const Real xend,Real& y,Real& h,const Real acc,
    const Real* pmt){
    const Real SAFETY=0.9;
    const Real PGROW=-0.2;
    const Real PSHRNK=-0.25;
    const Real ERRCON=1.89E-4;
    const Real TINY=1E-30;
    const int MAXSTP=30;
    const Real hmin=1E-12;
    Real ytemp,yerr;
    Real yscale;
    Real errmax;
    Real htemp;
    Real x=x0;
    for(int istp=0;istp<MAXSTP;istp++){
        yscale=fabs(y)+fabs(func(x,y,pmt)*h)+TINY;
        for(;;){
            RKCK(func,x,y,h,pmt,ytemp,yerr);
            errmax=0.0;
            errmax=max(errmax,fabs(yerr/yscale));
            errmax/=acc;
            if(errmax<=1.0){
                /*
                 * Succeed
                 */
                break;
            }
            htemp=SAFETY*h*pow(errmax,PSHRNK);
            h=(h>=0.0?max(htemp,0.1*h):min(htemp,0.1*h));
            if(h<hmin){
                cout<<"In RKQS, the stepsize underflow."<<endl;
                exit(1);
            }
        }
        // hdid=h;
        // if((x+h-xend)*(x+h-x0)>0.0){
        x+=h;
        y=ytemp;
        if(fabs(x-xend)<hmin){
            return;
        }
        h=(errmax>ERRCON?SAFETY*h*pow(errmax,PGROW):5.0*h);
        if(x+h>xend){
            h=xend-x;
        }
    }
} // RKQS

void RK45Vec(void func(const Real,const Real*,Real*,const Real*),
          Real& x,Real* YVec,const int n,const Real h,const Real* pmt) {
    Real k1[n],k2[n],k3[n],k4[n];
    Real YVecTmp[n];
    func(x,YVec,k1,pmt);
    for(int i=0;i<n;i++){
        YVecTmp[i]=YVec[i]+0.5*k1[i]*h;
    }
    func(x+0.5*h,YVecTmp,k2,pmt);
    for(int i=0;i<n;i++){
        YVecTmp[i]=YVec[i]+0.5*k2[i]*h;
    }
    func(x+0.5*h,YVecTmp,k3,pmt);
    for(int i=0;i<n;i++){
        YVecTmp[i]=YVec[i]+k3[i]*h;
    }
    func(x+h,YVecTmp,k4,pmt);
    for(int i=0;i<n;i++){
        YVec[i]+=h/6.0*(k1[i]+k2[i]*2.0+k3[i]*2.0+k4[i]);
    }
    x+=h;
} // RK45Vec

void RKAdaptiveVec(void func(const Real,const Real*,Real*,const Real*),
    const Real x0,const Real xend,Real* YVec,const int n,Real h,const Real acc,
    const Real* pmt) {
    // Real y0=y;
    Real x1,x2;
    Real YVec1[n],YVec2[n];
    Real x=x0;
    Real err,err_tmp;
    int count;
    const int MAXIT=100;
    while(fabs(x-xend)/xend>1E-9) {
    // while(x<xend) {
        // Take full step
        x1=x;
        for(int i=0;i<n;i++){
            YVec1[i]=YVec[i];
        }
        RK45Vec(func,x1,YVec1,n,h,pmt);
        // Take two half-steps
        x2=x;
        for(int i=0;i<n;i++){
            YVec2[i]=YVec[i];
        }
        RK45Vec(func,x2,YVec2,n,h/2.0,pmt);
        RK45Vec(func,x2,YVec2,n,h/2.0,pmt);
        // Use absolute truncation error
        // Take the max value of the array
        err=acc;
        for(int i=0;i<n;i++){
            err_tmp=fabs(YVec2[i]-YVec1[i])/15.0;
            err=err_tmp>err?err_tmp:err;
        }
        count=0;
        while(err>acc) {
            count++;
            if(count>MAXIT) {
                cout<<"MAX iteration is reached."<<endl;
                break;
            }
            // Calculate new step-length
            h*=pow(fabs(acc/err),0.2);
            h=(x+h)>xend?(xend-x):h;
            // Take full step
            x1=x;
            for(int i=0;i<n;i++){
                YVec1[i]=YVec[i];
            }
            RK45Vec(func,x1,YVec1,n,h,pmt);
            // Take two half-steps
            x2=x;
            for(int i=0;i<n;i++){
                YVec2[i]=YVec[i];
            }
            RK45Vec(func,x2,YVec2,n,h/2.0,pmt);
            RK45Vec(func,x2,YVec2,n,h/2.0,pmt);
            // Use absolute truncation error
            // Take the max value of the array
            err=acc;
            for(int i=0;i<n;i++){
                err_tmp=fabs(YVec2[i]-YVec1[i])/15.0;
                err=err_tmp>err?err_tmp:err;
            }
        }
        for(int i=0;i<n;i++){
            YVec[i]=YVec2[i];
        }
        x=x2;
        // cout<<x<<": ";
        // for(int i=0;i<n;i++){
            // cout<<YVec[i]<<",";
        // }
        // cout<<endl;
        // Calculate new step-length
        h*=pow(acc/err,0.2);
        // if(h<1E-12){
            // cout<<h<<endl;
            // cout<<acc<<endl;
            // cout<<err<<endl;
            // exit(1);
        // }
        h=(x+h)>xend?(xend-x):h;
    }
} // RKAdaptiveVec

#endif INTEGRATE_H