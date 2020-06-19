// fast and slow convolution routines, autoscaling and lamp shift
// Version 2006.03.26

#include "fsconv.h"
#include "fits2x.h"

/* rescaling -- old version. sum(fit)->sum(decay) */
void rescale(double *fit, double *decay, double *scale, int start, int stop)
{
    int i;
    double sumfit=0., sumcurve=0.;

    /* scaling */
    if (*scale==0.){
      for (i=start; i<=stop; i++){
        sumfit += fit[i];
        sumcurve += decay[i];
      }
      if (sumfit!=0.) *scale=sumcurve/sumfit;
     }
    for (i=start; i<=stop; i++)
      fit[i] *= *scale;
}

/* rescaling -- new version. scale = sum(fit*decay/w^2)/sum(fit^2/w^2) */
void rescale_w(double *fit, double *decay, double *w_sq, double *scale, int start, int stop)
{
    int i;
    double sumnom=0., sumdenom=0.;

    /* scaling */
    if (*scale==0.){
      for (i=start; i<=stop; i++){
        if (decay[i]!=0.) {
          sumnom += fit[i]*decay[i]/w_sq[i];
          sumdenom += fit[i]*fit[i]/w_sq[i];
        }
      }
      if (sumdenom!=0.) *scale = sumnom/sumdenom;
     }
    for (i=start; i<=stop; i++)
      fit[i] *= *scale;

}

/* rescaling -- new version + background. scale = sum(fit*decay/w^2)/sum(fit^2/w^2) */
void rescale_w_bg(double *fit, double *decay, double *w_sq, double bg, double *scale, int start, int stop)
{
    int i;
    double sumnom=0., sumdenom=0.;

    /* scaling */
    if (*scale==0.){
      for (i=start; i<=stop; i++){
        if (decay[i]!=0.) {
          sumnom += fit[i]*(decay[i]-bg)/w_sq[i];
          sumdenom += fit[i]*fit[i]/w_sq[i];
        }
      }
      if (sumdenom!=0.) *scale = sumnom/sumdenom;
     }
    for (i=start; i<=stop; i++)
      fit[i] *= *scale;

}

/* fast convolution */
void fconv(double *fit, double *x, double *lamp, int numexp, int start, int stop)
{
    int ne,i;
    double fitcurr, expcurr, deltathalf = DELTA_T*0.5;

    for (i=0;i<=stop;i++) fit[i]=0; 

    /* convolution */
    for (ne=0; ne<numexp; ne++) {
      expcurr = exp(-DELTA_T/x[2*ne+1]);
      fitcurr = 0;
      for (i=1; i<=stop; i++){
        fitcurr=(fitcurr + deltathalf*lamp[i-1])*expcurr + deltathalf*lamp[i];
        fit[i] += fitcurr*x[2*ne];
      }
    }

}

/* fast convolution, high repetition rate */
void fconv_per(double *fit, double *x, double *lamp, int numexp, int start, int stop,
           int n_points, double period)
{
    int ne, i, lamp_start = 0, 
      stop1, period_n = (int)ceil(period/DELTA_T-0.5);
    double fitcurr, expcurr, tail_a, deltathalf = DELTA_T*0.5;

    while (lamp[lamp_start++]==0);
    for (i=0;i<=stop;i++) fit[i]=0;

    stop1 = (period_n+lamp_start > n_points-1) ? n_points-1 : period_n+lamp_start;

    /* convolution */
    for (ne=0; ne<numexp; ne++) {
      expcurr = exp(-DELTA_T/x[2*ne+1]);
      tail_a = 1./(1.-exp(-period/x[2*ne+1]));
      fitcurr = 0;
      for (i=1; i<=stop1; i++){
        fitcurr=(fitcurr + deltathalf*lamp[i-1])*expcurr + deltathalf*lamp[i];
        fit[i] += fitcurr*x[2*ne];
      }
      fitcurr *= exp(-(period_n - stop1 + start)*DELTA_T/x[2*ne+1]);
      for (i=start; i<=stop; i++){
        fitcurr *= expcurr;
        fit[i] += fitcurr*x[2*ne]*tail_a;
      }
    }

}

/* fast convolution, high repetition rate, with convolution stop for Paris */
void fconv_per_cs(double *fit, double *x, double *lamp, int numexp, int stop,
           int n_points, double period, int conv_stop)
{
    int ne, i, 
      stop1, period_n = (int)ceil(period/DELTA_T-0.5);
    double fitcurr, expcurr, tail_a, deltathalf = DELTA_T*0.5;

    for (i=0; i<=stop; i++) fit[i]=0;
    stop1 = (period_n > n_points-1) ? n_points-1 : period_n;

    /* convolution */
    for (ne=0; ne<numexp; ne++) {
      expcurr = exp(-DELTA_T/x[2*ne+1]);
      tail_a = 1./(1.-exp(-period/x[2*ne+1]));
      fitcurr = 0.;
      fit[0] += deltathalf*lamp[0]*(expcurr + 1.)*x[2*ne];
      for (i=1; i<=conv_stop; i++) {
        fitcurr=(fitcurr + deltathalf*lamp[i-1])*expcurr + deltathalf*lamp[i];
        fit[i] += fitcurr*x[2*ne];
      }
      for (; i<=stop1; i++) {
        fitcurr=fitcurr*expcurr;
        fit[i] += fitcurr*x[2*ne];
      }
      fitcurr *= exp(-(period_n - stop1)*DELTA_T/x[2*ne+1]);
      for (i=0; i<=stop; i++) {
        fitcurr *= expcurr;
        fit[i] += fitcurr*x[2*ne]*tail_a;
      }
    }
}

/* fast convolution with reference compound decay */
void fconv_ref(double *fit, double *x, double *lamp, int numexp, int start, int stop, double tauref)
{
    int ne,i;
    double fitcurr, expcurr, deltathalf = DELTA_T*0.5, correct_a, sum_a=0;

    for (i=0;i<=stop;i++) fit[i]=0; 

    /* convolution */
    for (ne=0; ne<numexp; ne++) {
      expcurr = exp(-DELTA_T/x[2*ne+1]);
      correct_a = x[2*ne]*(1/tauref-1/x[2*ne+1]);
      sum_a += x[2*ne];
      fitcurr = 0;
      for (i=1; i<=stop; i++){
        fitcurr=(fitcurr + deltathalf*lamp[i-1])*expcurr + deltathalf*lamp[i];
        fit[i] += fitcurr * correct_a;
      }
    }
    for (i=1; i<=stop; i++) fit[i] += lamp[i] * sum_a;

}

/* slow convolution */
void sconv(double *fit, double *p, double *lamp, int start, int stop)
{
    int i,j;

    /* convolution */
    for (i=start; i<=stop; i++){
       fit[i] = 0.5 * lamp[0] * p[i];
       for (j=1; j<i; j++) fit[i] += lamp[j] * p[i-j];
       fit[i] += 0.5 * lamp[i] * p[0];
       fit[i]=fit[i]*DELTA_T;
    }
    fit[0] = 0;

}


/* shifting lamp */
void shift_lamp(double *lampsh, double *lamp, double ts, int n_points)
{
    int tsint = (int)(floor(ts));
    double tsdbl = ts-(double)tsint;
    int out_left=0, out_right=0, j; 

    if (tsint<0) out_left = -tsint;
    if (tsint+1>0) out_right = tsint+1;

    for(j=0; j<out_left; j++) lampsh[j]=0;
    for(j=out_left; j<(n_points-out_right); j++) 
       lampsh[j]=lamp[j+tsint]*(1-tsdbl)+lamp[j+tsint+1]*(tsdbl);
    for(j=(n_points-out_right); j<n_points; j++) lampsh[j]=0;

}

