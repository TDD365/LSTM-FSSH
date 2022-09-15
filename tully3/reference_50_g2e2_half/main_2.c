#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#define EXC
#define max(a,b) ((a)>(b)?(a):(b))
//lawuwuwuwuwu

int diamodel(double x, double ener[3], double grad[3])
{
    double y;
    double a=6.0, b=0.1, c=0.9, a0=10000.0;

    ener[0] = a / a0;
    ener[1] = -1.0 * ener[0];
    grad[0] = grad[1] = 0.0;

    if (x > 0.0) {
        y = exp(-1.0 * c * x);
        ener[2] = b * (2.0 - y);
        grad[2] = b * c * y;
    }
    else {
        y = exp(c * x);
        ener[2] = b * y;
        grad[2] = b * c * y;
    }

    return 0;
}

int diamodel2(double x, double ener[3], double grad[3])
{
    double y;
    double a=0.1, b=0.28, c=0.015, d=0.06, e=0.05;

    ener[0] = 0.0;
    ener[1] = - a * exp(-1.0 * b * x * x) + e;
    ener[2] =   c * exp(-1.0 * d * x * x);
    
    grad[0] = 0.0;
    grad[1] =   2.0 * a * b * x * exp(-1.0 * b * x * x);
    grad[2] = - 2.0 * c * d * x * exp(-1.0 * d * x * x);

    return 0;
}

int diamodel1(double x, double ener[3], double grad[3])
{
    double y;
    double a=0.01, b=1.6, c=0.005, d=1.0;
    
    if (x > 0.0) {
        y = exp(-1.0 * b * x);
        ener[0] = a * (1.0 - y);
        grad[0] = a * b * y;
        ener[1] = -1.0 * ener[0];
        grad[1] = -1.0 * grad[0];
    }
    else {
        y = exp(b * x);
        ener[0] = -1.0 * a * (1.0 - y);
        grad[0] = a * b * y;
        ener[1] = -1.0 * ener[0];
        grad[1] = -1.0 * grad[0];
    }

    ener[2] = c * exp(-1.0 * d * x * x);
    grad[2] = - 2.0 * c * d * x * exp(-1.0 * d * x * x);

    return 0;
}

int adener(double xt, double ener[2], double grad[3])
{
    double x, y, z, dhalm[3], ghalm[3];

    diamodel(xt, dhalm, ghalm);

    x = (dhalm[0] - dhalm[1]) * (dhalm[0] - dhalm[1]) + 4.0 * dhalm[2] * dhalm[2];
    y = sqrt(x);
    z = (dhalm[0] - dhalm[1]) * (ghalm[0] - ghalm[1]) + 4.0 * dhalm[2] * ghalm[2];
    ener[0] = 0.5 * (dhalm[0] + dhalm[1] - y);
    ener[1] = 0.5 * (dhalm[0] + dhalm[1] + y);
    grad[0] = 0.5 * (ghalm[0] + ghalm[1] - z / y);
    grad[1] = 0.5 * (ghalm[0] + ghalm[1] + z / y);
    grad[2] = (ghalm[2] * (dhalm[0] - dhalm[1]) - dhalm[2] * (ghalm[0] - ghalm[1])) / x;

    return 0;
}

void Func(double v,double complex rho[],double complex d[],double g,double e0,double e1)
{
    d[0*2+0] = -v*(g*(rho[0*2+1]+rho[1*2+0]));
    d[0*2+1] = v*(g*(rho[0*2+0]-rho[1*2+1]))-I*rho[0*2+1]*(e0-e1);
    d[1*2+0] = v*(g*(rho[0*2+0]-rho[1*2+1]))-I*rho[1*2+0]*(e1-e0);
    d[1*2+1] = v*(g*(rho[0*2+1]+rho[1*2+0]));
    return ;
}

void RKT(double v,double complex rho[],double h,double g,double gt,double e0,double et0,double e1,double et1)
{
    extern void Func();
    int i,j,m;
    double a[4],gm[3],em[3],eem[3];
    double complex b[4],d[4],z[4];

    a[0]=h/2.0;
    a[1]=h/2.0;
    a[2]=h;
    a[3]=h;

    gm[0]=1.0/2.0*(g+gt);
    gm[1]=1.0/2.0*(g+gt);
    gm[2]=gt;

    em[0]=1.0/2.0*(e0+et0);
    em[1]=1.0/2.0*(e0+et0);
    em[2]=et0;

    eem[0]=1.0/2.0*(e1+et1);
    eem[1]=1.0/2.0*(e1+et1);
    eem[2]=et1;

    Func(v,rho,d,g,e0,e1);
    for (i=0; i<=1; i++) {
        for(j=0; j<=1;j++) {
            b[i*2+j]=rho[i*2+j];
            z[i*2+j]=rho[i*2+j];
        }
    }
    for (m=0; m<=2; m++) {
        for (i=0; i<=1; i++) {
            for(j=0; j<=1;j++) {
                rho[i*2+j]=z[i*2+j]+a[m]*d[i*2+j];
                b[i*2+j]=b[i*2+j]+a[m+1]*d[i*2+j]/3.0;
            }
        }
        g = gm[m];
        e0 = em[m];
        e1 = eem[m];
        Func(v,rho,d,g,e0,e1);
    }
    for(i=0; i<=1; i++) {
        for(j=0; j<=1;j++)
            rho[i*2+j]=b[i*2+j]+h*d[i*2+j]/6.0;
    }
    return ;
}

int main()
{
    int istep, istate, istate_x, i, j;
    int T1, T2, R1, R2, inv;
    double ip, p, t, dt, m, x, xt, v, vt, a, at, g, gt, e0, et0, e1, et1;
    double random, P_ij, g_ij, a_ij, b_ij, c_ij;
    double g_fix,e0_fix,e1_fix,a_fix;
    double ener[2], grad[3];
    double complex rho[4];

    double P_T1,P_T2,P_R1,P_R2;

    FILE *fp,*fpt;
    fp = fopen("result2.txt","w");
    if((fp = fopen("result2.txt","w"))==NULL) {
        printf("can not open the file!\n");
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "%s","p   T1     T2     R1     R2 \n");
    fpt = fopen("result_traj_2.txt","w");
    
    srand((unsigned)time( NULL ));
    for(ip = 3.0; ip <= 32.0; ip=ip+1.0 ) {
    //for(ip = -3.5; ip <= 0.2; ip=ip+0.1 ) {
        
        dt = 1.0;
        m = 2000.0;
        x = -20.0;
        adener(x, ener, grad);
        g = grad[2];
        e0 = ener[0];
        e1 = ener[1];
        a = -grad[0]/m;
        g_fix = g;
        e0_fix = e0;
        e1_fix = e1;
        a_fix = a;

        //p = sqrt(2.0*m*exp(ip));
        p = ip;
        fprintf(fpt,"VELOC %6.4f\n\n", ip);
        
        T1 = 0;
        T2 = 0;
        R1 = 0;
        R2 = 0;
        inv = 0;   //invalid

        for(i=0;i<50;i++) {
            
            fprintf(fpt,"TRAJ %d\n", i);
            
            t = 0.0;
            istate = 0;
            istep = 0;
            rho[0*2+0] = 1.0;
            rho[0*2+1] = 0.0;
            rho[1*2+0] = 0.0;
            rho[1*2+1] = 0.0;
            x = -20.0;
            v = (double)(p)/m;
            g = g_fix;
            e0 = e0_fix;
            e1 = e1_fix;
            a = a_fix;

            while(x<0 && v>0) {
                
                istep++;
                xt = x+v*dt+0.5*a*dt*dt;
                adener(xt, ener, grad);
                gt = grad[2];
                et0 = ener[0];
                et1 = ener[1];
                at = -grad[istate]/m;
                vt = v+0.5*dt*(a+at);
                RKT(v,rho,dt,g,gt,e0,et0,e1,et1);
                istate_x = istate;

                fprintf(fpt,"%12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %d %12.9f %12.9f ", x, v, e0, e1, creal(rho[0]), creal(rho[1]), cimag(rho[1]), istate, fabs(v)*grad[0]/e0/e0/15.0, fabs(v)*grad[1]/e1/e1/15.0);

                for(j=0;j<=1;j++) {
                    if(j!=istate_x) {
                        if(j<istate_x)
                            grad[2] = -grad[2];
                        random = ((double)rand()/RAND_MAX);
                        P_ij = 2.0*dt*creal(rho[istate*2+j]*vt*grad[2])/rho[istate*2+istate];
                        g_ij = max(P_ij,0.0);
                        if( random < g_ij ) {
                            a_ij = 0.5*grad[2]*grad[2]/m;
                            b_ij = vt*grad[2];
                            if( b_ij*b_ij+4.0*a_ij*(ener[istate]-ener[j]) >= 0 ) {
                                if(b_ij < 0)
                                    c_ij = (b_ij+sqrt(b_ij*b_ij+4.0*a_ij*(ener[istate]-ener[j])))/(2.0*a_ij);  // a scaling factor to be determined
                                else
                                    c_ij = (b_ij-sqrt(b_ij*b_ij+4.0*a_ij*(ener[istate]-ener[j])))/(2.0*a_ij);
                                vt = vt-c_ij*grad[2]/m;
                                at = -grad[j]/m;
                                istate = j;
                            }
                        }
                    }
                }

                fprintf(fpt,"%12.9f \n", g_ij);

                v = vt;
                x = xt;
                a = at;
                t = t+dt;
                g = gt;
                e0 = et0;
                e1 = et1;
            }
        if ( istep<1000000 && x>25 && istate==0 )
            T1++;
        else if ( istep<1000000 && x>25 && istate==1 )
            T2++;
        else if ( istep<1000000 && x<-25 && istate==0 )
            R1++;
        else if ( istep<1000000 && x<-25 && istate==1 )
            R2++;
        else
            inv++;
        }
        P_T1 = (double)(T1)/(double)(50-inv);
        P_T2 = (double)(T2)/(double)(50-inv);
        P_R1 = (double)(R1)/(double)(50-inv);
        P_R2 = (double)(R2)/(double)(50-inv);
        fprintf(fp,"%6.4f %6.4f %6.4f %6.4f %6.4f \n", ip, P_T1, P_T2, P_R1, P_R2);
        fflush(fp);
        fflush(fpt);

    }

    fclose(fp);
    fclose(fpt);

}


