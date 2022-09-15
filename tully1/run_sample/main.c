#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <time.h>
#define EXC
#define max(a,b) ((a)>(b)?(a):(b))
//lawuwuwuwuwu
#define LSTMON

// LSTM
#include "./include/k2c_include.h"
#include "./include/k2c_tensor_include.h"
#include "LSTM_test.h"

/* MODELS */
int diamodel3(float x, float ener[3], float grad[3])
{
    float y;
    float a=6.0, b=0.1, c=0.9, a0=10000.0;

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

int diamodel2(float x, float ener[3], float grad[3])
{
    float a=0.1, b=0.28, c=0.015, d=0.06, e=0.05;

    ener[0] = 0.0;
    ener[1] = - a * exp(-1.0 * b * x * x) + e;
    ener[2] =   c * exp(-1.0 * d * x * x);

    grad[0] = 0.0;
    grad[1] =   2.0 * a * b * x * exp(-1.0 * b * x * x);
    grad[2] = - 2.0 * c * d * x * exp(-1.0 * d * x * x);

    return 0;
}

int diamodel(float x, float ener[3], float grad[3])
{
    float y;
    float a=0.01, b=1.6, c=0.005, d=1.0;

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

int adener(float xt, float ener[2], float grad[3])
{
    float x, y, z, dhalm[3], ghalm[3];

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

void Func(float v,float complex rho[],float complex d[],float g,float e0,float e1)
{
    d[0*2+0] = -v*(g*(rho[0*2+1]+rho[1*2+0]));
    d[0*2+1] = v*(g*(rho[0*2+0]-rho[1*2+1]))-I*rho[0*2+1]*(e0-e1);
    d[1*2+0] = v*(g*(rho[0*2+0]-rho[1*2+1]))-I*rho[1*2+0]*(e1-e0);
    d[1*2+1] = v*(g*(rho[0*2+1]+rho[1*2+0]));
    return ;
}

/* RK4 Process */
void RKT(float v,float complex rho[],float h,float g,float gt,float e0,float et0,float e1,float et1)
{
    int i,j,m;
    float a[4],gm[3],em[3],eem[3];
    float complex b[4],d[4],z[4];

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

float * get_array(k2c_tensor *LSTM_output) {

    return LSTM_output->array;

}

/* Calling LSTM */
float * call_LSTM(float LSTM_input[700]) {

    // First, make tensors
    float LSTM_output_fake[3] = {0};
    k2c_tensor lstm_input_work           = {&LSTM_input[0],2,700,{100,  7,  1,  1,  1}};
    k2c_tensor lstm_output_work          = {&LSTM_output_fake[0],1,3,{3,1,1,1,1}};

    LSTM_test(&lstm_input_work,&lstm_output_work);

    return get_array(&lstm_output_work);

}

/* Main Process */
int main()
{
    int nrun=16, istep, istate, istate_x, i, j;
    int T1, T2, R1, R2, inv;
    float ip, p, t, dt, m, x, xt, v, vt, a, at, g, gt, e0, et0, e1, et1;
    float random, P_ij, g_ij, a_ij, b_ij, c_ij;
    float g_fix,e0_fix,e1_fix,a_fix;
    float ener[2], grad[3];
    float complex rho[4];
#ifdef LSTMON
    int k ,l;
    float rho_old[4];
#endif

    float P_T1,P_T2,P_R1,P_R2;

    /*
    This is a matrix that LSTM needs, including

    X,V,E0,E1,RERHO0,RERHO1,IMRHO1(,ACTIVE,POSS)

    */
    int LSTM_count=0;
    int LSTM_LEARN=100;
    int LSTM_IN=7;
    int LSTM_OUT=3;
    float LSTM_input[LSTM_LEARN*LSTM_IN];
    float LSTM_output[LSTM_LEARN*LSTM_OUT];

    // File for trajectory & collective result
    FILE *fp, *fpt;
    fp  = fopen("result.txt","w");
    fpt = fopen("result_traj.txt","w");
    // Table Title
    fprintf(fp, "%s","p   T1     T2     R1     R2 \n");

    // Random Seed
    srand((unsigned)time(NULL));

    // Go!
    printf("Tully's Model in Surface Hopping, LSTM Version, 8/27/2021\n");
    printf("C Programming by JiaLY\n");
    printf("C++/LSTM compatibility done by TangDD\n\n");

    for(ip = 3.0; ip <= 32.0; ip=ip+1.0 ) {

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
        p=ip;
        printf("VELOC %6.4f\n\n", ip);
        //fprintf(fpt,"VELOC %6.4f\n\n", ip);

        T1 = 0;
        T2 = 0;
        R1 = 0;
        R2 = 0;
        inv = 0;   //invalid

        for(i=0;i<nrun;i++) {

            //printf("TRAJ %d\n", i);

            // Initialize LSTM Matrix
            LSTM_count=0;
            memset(LSTM_input, 0, sizeof(LSTM_input));
            memset(LSTM_output, 0, sizeof(LSTM_output));

            t = 0.0;
            istate = 0;
            istep = 0;
            rho[0*2+0] = 1.0;
            rho[0*2+1] = 0.0;
            rho[1*2+0] = 0.0;
            rho[1*2+1] = 0.0;
            x = -20.0;
            v = (float)(p)/m;
            g = g_fix;
            e0 = e0_fix;
            e1 = e1_fix;
            a = a_fix;

            while(fabs(x)<=25 && istep<1000000) {

#ifdef LSTMON
                // LSTM: We will record at step start:
                // X,V,E0,E1,RERHO0,RERHO1,IMRHO1
                if (LSTM_count<LSTM_LEARN) {
                    LSTM_input[LSTM_count*LSTM_IN+0] = x/20.0;
                    LSTM_input[LSTM_count*LSTM_IN+1] = (v*2000-17.0)/30.0;
                    LSTM_input[LSTM_count*LSTM_IN+2] = e0*50.0;
                    LSTM_input[LSTM_count*LSTM_IN+3] = e1*50.0;
                    LSTM_input[LSTM_count*LSTM_IN+4] = creal(rho[0])-0.5;
                    LSTM_input[LSTM_count*LSTM_IN+5] = creal(rho[1]);
                    LSTM_input[LSTM_count*LSTM_IN+6] = cimag(rho[1]);
                } else {
                    for (k=1;k<LSTM_LEARN;k++) {
                        for (l=0;l<LSTM_IN;l++) {
                            LSTM_input[(k-1)*LSTM_IN+l] = LSTM_input[k*LSTM_IN+l];
                        }
                    }
                    // Now overwrite!
                    LSTM_input[(LSTM_LEARN-1)*LSTM_IN+0] = x/20.0;
                    LSTM_input[(LSTM_LEARN-1)*LSTM_IN+1] = (v*2000-17.0)/30.0;
                    LSTM_input[(LSTM_LEARN-1)*LSTM_IN+2] = e0*50.0;
                    LSTM_input[(LSTM_LEARN-1)*LSTM_IN+3] = e1*50.0;
                    LSTM_input[(LSTM_LEARN-1)*LSTM_IN+4] = creal(rho[0])-0.5;
                    LSTM_input[(LSTM_LEARN-1)*LSTM_IN+5] = creal(rho[1]);
                    LSTM_input[(LSTM_LEARN-1)*LSTM_IN+6] = cimag(rho[1]);
                }
                LSTM_count += 1;
#endif
                istep++;
                xt = x+v*dt+0.5*a*dt*dt;
                adener(xt, ener, grad);
                gt = grad[2];
                et0 = ener[0];
                et1 = ener[1];
                at = -grad[istate]/m;
                vt = v+0.5*dt*(a+at);

                rho_old[0*2+0] = rho[0*2+0];
                rho_old[0*2+1] = rho[0*2+1];
                rho_old[1*2+0] = rho[1*2+0];
                rho_old[1*2+1] = rho[1*2+1];

#ifdef LSTMON
                // RK4 Process, which is to be replaced by LSTM integrator.
                if (LSTM_count<LSTM_LEARN) {
                    // This is Traditional RK4.
                    RKT(v,rho,dt,g,gt,e0,et0,e1,et1);
                } else {
                    // Calling LSTM.
                    memcpy(LSTM_output,call_LSTM(LSTM_input),sizeof(LSTM_output));
                    rho[0]= (0.5+LSTM_output[0]) + 0.0*I;
                    rho[1]=      LSTM_output[1]  + LSTM_output[2]*I;
                    rho[2]=      LSTM_output[1]  - LSTM_output[2]*I;
                    rho[3]= (0.5-LSTM_output[0]) + 0.0*I;
                }
#else
                RKT(v,rho,dt,g,gt,e0,et0,e1,et1);
#endif
                istate_x = istate;

                //if (istep%100==0) fprintf(fpt,"%12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %d ", x, v, e0, e1, creal(rho[0]), creal(rho[1]), cimag(rho[1]), istate);
                //if (istep%1==0)
                //printf("%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %d \n", x, v, e0, e1, creal(rho[0]), creal(rho[1]), cimag(rho[1]), istate);

                for(j=0;j<=1;j++) {
                    if(j!=istate_x) {
                        if(j<istate_x)
                            grad[2] = -grad[2];
                        random = ((float)rand()/RAND_MAX);
#ifdef LSTMON
                        //P_ij = 2.0*dt*creal(rho[istate*2+j]*vt*grad[2])/rho[istate*2+istate];
                        if (LSTM_count<LSTM_LEARN) {
                            // This is Traditional probability.
                            P_ij = 2.0*dt*creal(rho[istate*2+j]*vt*grad[2])/rho[istate*2+istate];
                        } else {
                            // In LSTM, we only have numeral gradient of density
                            //P_ij = 2.0*(creal(rho_old[3*istate])-creal(rho[3*istate]))/(creal(rho_old[3*istate])+creal(rho[3*istate]));
                            P_ij = (creal(rho_old[3*istate])-creal(rho[3*istate]))/creal(rho_old[3*istate]);
                        }

#else
                        P_ij = 2.0*dt*creal(rho[istate*2+j]*vt*grad[2])/creal(rho[istate*2+istate]);
#endif
                        g_ij = max(P_ij,0.0);
#ifdef LSTMON
                        // TangDD, modified E+
                        if( random < g_ij ) {
                            // The kinetic energy;
                            a_ij = 0.5 * vt * vt * m;
                            // Free potential energy, j is target state
                            b_ij = ener[j] - ener[istate];
                            // Check if energy can converge?
                            if( a_ij - b_ij >= 0 ) {
                                // scaling factor
                                c_ij = sqrt(1.0 - b_ij / a_ij);
                                vt = vt * c_ij;
                                istate = j;
                            }
                        }
#else
                        // Standard E+
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
#endif
                    }
                }

                //if (istep%100==0) fprintf(fpt,"%12.9f \n", g_ij);

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
        P_T1 = (float)(T1)/(float)(nrun-inv);
        P_T2 = (float)(T2)/(float)(nrun-inv);
        P_R1 = (float)(R1)/(float)(nrun-inv);
        P_R2 = (float)(R2)/(float)(nrun-inv);
        fprintf(fp,"%6.4f %6.4f %6.4f %6.4f %6.4f \n", ip, P_T1, P_T2, P_R1, P_R2);
        fflush(fp);
        fflush(fpt);

    }

    fclose(fp);
    fclose(fpt);

    return 0;

}


