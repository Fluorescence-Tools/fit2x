#include <iomanip>      // std::setprecision
#include "fits2x.h"
#include "fsconv.h"
#include "statistics.h"

// normalization
static double Sp, Ss, Bp, Bs;
static double Bexpected;    // expected <B> corresponding to the mean Bg signal
static int fixedrho = 0;
static int softbifl = 0;
static int p2s_twoIstar = 0;
static int firstcall = 1;
static double penalty = 0.;



void compute_signal_and_background(MParam *p) {
#if VERBOSE_FIT2X
    std::cout << "COMPUTE_SIGNAL_AND_BACKGROUND" << std::endl;
#endif
    // total signal and background
    int i;
    LVI32Array *expdata = *(p->expdata);
    int Nchannels = expdata->length / 2;
    LVDoubleArray *irf = *(p->irf), *bg = *(p->bg),
            *corrections = *(p->corrections), *M = *(p->M);

    // total signal and background
    Sp = 0.;
    Ss = 0.;
    Bp = 0.;
    Bs = 0.;
    for (i = 0; i < Nchannels; i++) {
        Sp += expdata->data[i];
        Bp += bg->data[i];
    }
    for (; i < 2 * Nchannels; i++) {
        Ss += expdata->data[i];
        Bs += bg->data[i];
    }
    double B = std::max(1.0, Bp + Bs);
    Bp *= (Sp + Ss) / B;
    Bs *= (Sp + Ss) / B;
#if VERBOSE_FIT2X
    std::cout << "-- Bp, Bs: " << Bp << ", " << Bs << std::endl;
    std::cout << "-- Sp, Ss: " << Sp << ", " << Ss << std::endl;
#endif
}


void normM(double *M, int Nchannels) {
    int i;
    double s = 0., Sexp = Sp + Ss;
#if VERBOSE_FIT2X
    std::cout << "NORMALIZE: normM" << std::endl;
    std::cout << "-- Nchannels: " << Nchannels << std::endl;
    std::cout << "-- Sp: " << Sp << std::endl;
    std::cout << "-- Ss: " << Ss << std::endl;
    std::cout << "-- Sexp: " << Ss << std::endl;
    std::cout << "-- s: " << s << std::endl;
#endif
    for (i = 0; i < 2 * Nchannels; i++) s += M[i];
    for (i = 0; i < 2 * Nchannels; i++) M[i] *= Sexp / s;

}

// if already normalized to s:
void normM(double *M, double s, int Nchannels) {
    int i;
    double Sexp = Sp + Ss;
#if VERBOSE_FIT2X
    std::cout << "NORMALIZE: normM - already normalized to s" << std::endl;
    std::cout << "-- Nchannels: " << Nchannels << std::endl;
    std::cout << "-- s: " << s << std::endl;
    std::cout << "-- Sp: " << Sp << std::endl;
    std::cout << "-- Ss: " << Ss << std::endl;
    std::cout << "-- Sexp: " << Ss << std::endl;
#endif
    for (i = 0; i < 2 * Nchannels; i++) M[i] *= Sexp / s;
}

void normM_p2s(double *M, int Nchannels) {
    int i;
    double s = 0.;

    for (i = 0; i < Nchannels; i++) s += M[i];
    for (i = 0; i < Nchannels; i++) M[i] *= Sp / s;

    s = 0.;
    for (i = Nchannels; i < 2 * Nchannels; i++) s += M[i];
    for (i = Nchannels; i < 2 * Nchannels; i++) M[i] *= Ss / s;
}


//////////////////////////////////////////// FIT23 ////////////////////////////////////////////


void correct_input23(double *x, double *xm, LVDoubleArray *corrections, int return_r) {
    double r, g, Fp, Fs, l1, l2;
#if VERBOSE_FIT2X
    std::cout << "CORRECT_INPUT23" << std::endl;
#endif
    xm[0] = x[0];
    if (xm[0] < 0.001) {
        xm[0] = 0.001;    // tau > 0
        penalty = -x[0];
    } else penalty = 0.;

    if (x[1] < 0.) xm[1] = 0.;        // 0 < gamma < 0.999
    else if (x[1] > 0.999) xm[1] = 0.999;
    else xm[1] = x[1];

    xm[2] = x[2];
    // anisotropy
    g = corrections->data[1];
    l1 = corrections->data[2];
    l2 = corrections->data[3];
#if VERBOSE_FIT2X
    std::cout << "-- Correction factors:" << std::endl;
    std::cout << "-- g-factor: " << g << std::endl;
    std::cout << "-- l1, l2: " << l1 << ", " << l2 << std::endl;
#endif
    Fp = (Sp - xm[1] * Bp) / (1. - xm[1]);
    Fs = (Ss - xm[1] * Bs) / (1. - xm[1]);
    r = (Fp - g * Fs) / (Fp * (1. - 3. * l2) + (2. - 3. * l1) * g * Fs);
#if VERBOSE_FIT2X
    std::cout << "-- Signals: " << std::endl;
    std::cout << "-- Bp, Bs: " << Bp << ", " << Bs << std::endl;
    std::cout << "-- Sp, Ss: " << Sp << ", " << Ss << std::endl;
    std::cout << "-- Fp, Fs: " << Fp << ", " << Fs << std::endl;
    std::cout << "-- r: " << r << std::endl;
#endif
    if (!fixedrho) {
        xm[3] = x[0] / (x[2] / r - 1.);        // rho = tau/(r0/r-1)
        if (xm[3] < 1.e-4) xm[3] = 1.e-4;    // rho > 0
        x[3] = xm[3];
    } else xm[3] = x[3];

    if (return_r) {
        x[7] = (Sp - g * Ss) / (Sp * (1. - 3. * l2) + (2. - 3. * l1) * g * Ss);
        x[6] = r;
    }
#if VERBOSE_FIT2X
    std::cout << "-- Corrected parameters:" << std::endl;
    std::cout << "-- tau / tau_c: " << x[0] << "/" << xm[0] << std::endl;
    std::cout << "-- gamma / gamma_c: " << x[1] << "/" << xm[1] << std::endl;
    std::cout << "-- r0 / r0_c: " << x[2] << "/" << xm[2] << std::endl;
    std::cout << "-- rho / rho_c: " << x[3] << "/" << xm[3] << std::endl;
#endif
}


double targetf23(double *x, void *pv) {
    double w, xm[4], Bgamma;
    MParam *p = (MParam *) pv;
    LVI32Array *expdata = *(p->expdata);
    int Nchannels = expdata->length / 2;
    LVDoubleArray
        *irf = *(p->irf),
        *bg = *(p->bg),
        *corrections = *(p->corrections),
        *M = *(p->M);
#if VERBOSE_FIT2X
    std::cout << "COMPUTING TARGET23" << std::endl;
    std::cout << "-- Nchannels: " << Nchannels << std::endl;
#endif
    correct_input23(x, xm, corrections, 0);
    modelf23(xm, irf->data, bg->data, Nchannels, p->dt, corrections->data, M->data);
    normM(M->data, 1., Nchannels);
    if (p2s_twoIstar) w = Wcm_p2s(expdata->data, M->data, Nchannels);
    else w = Wcm(expdata->data, M->data, Nchannels);
    if (softbifl && (Bexpected > 0.)) {
        Bgamma = xm[1] * (Sp + Ss);
        w -= Bgamma * log(Bexpected) - loggammaf(Bgamma + 1.);
    }
    double v = w / Nchannels + penalty;
//#if VERBOSE_FIT2X
//    std::cout << "xm:" ;
//    for(int i=0; i<8;i++) std::cout << xm[i] << " ";
//    std::cout << std::endl;
//
//    std::cout << "corrections:" ;
//    for(int i=0; i<corrections->length;i++) std::cout << corrections->data[i] << " ";
//    std::cout << std::endl;
//
//    std::cout << "irf:" ;
//    for(int i=0; i<irf->length;i++) std::cout << irf->data[i] << " ";
//    std::cout << std::endl;
//
//    std::cout << "bg:" ;
//    for(int i=0; i<bg->length;i++) std::cout << bg->data[i] << " ";
//    std::cout << std::endl;
//
//    std::cout << "Data:" ;
//    for(int i=0; i<expdata->length;i++) std::cout << expdata->data[i] << " ";
//    std::cout << std::endl;
//
//    std::cout << "Model:" ;
//    for(int i=0; i<M->length;i++) std::cout << std::setprecision(2) << M->data[i] << " ";
//    std::cout << std::endl;
//
//    std::cout << "g" << std::endl;
//    std::cout << v << std::endl;
//#endif
    return v;
}


double fit23(double *x, short *fixed, MParam *p) {
#if VERBOSE_FIT2X
    std::cout << "FIT23" << std::endl;
    std::cout << "-- Initial parameters / fixed: " << std::endl;
    std::cout << "-- tau: " << x[0] << " / " << fixed[0] << std::endl;
    std::cout << "-- gamma: " << x[1] << " / " << fixed[1] << std::endl;
    std::cout << "-- r0: " << x[2] << " / " << fixed[2] << std::endl;
    std::cout << "-- rho: " << x[3] << " / " << fixed[3] << std::endl;
    std::cout << "-- Soft BIFL scatter fit?: " << x[4] << std::endl;
    std::cout << "-- 2I*: P+2S?: " << x[5] << std::endl;
    std::cout << "-- r Scatter (output only): " << x[6] << std::endl;
    std::cout << "-- r Experimental (output only): " << x[7] << std::endl;
#endif
    double tIstar, xm[4];
    int info = -1;

    if (firstcall) init_fact();
    firstcall = 0;

    softbifl = (x[4] < 0.);
    p2s_twoIstar = (x[5] > 0.);

    LVI32Array *expdata = *(p->expdata);
    int Nchannels = expdata->length / 2;
    LVDoubleArray *corrections = *(p->corrections), *M = *(p->M);

    compute_signal_and_background(p);
    Bexpected = x[1] * (Sp + Ss);
    fixedrho = fixed[3];

    bfgs bfgs_o(targetf23, 4);

    bfgs_o.fix(1);    // gamma
    bfgs_o.fix(2);    // r0
    bfgs_o.fix(3);    // rho is set in targetf23

    // pre-fit with fixed gamma
    //  bfgs_o.maxiter = 20;
    if(!fixed[0]){
#if VERBOSE_FIT2X
        std::cout << "-- pre-fitting..." << std::endl;
#endif
        bfgs_o.minimize(x, p);
    }else {
        bfgs_o.fix(0);
    }

    // fit with free gamma
    // bfgs_o.maxiter = 100;
    if (!fixed[1] && (x[4] <= 0.)) {
        bfgs_o.free(1);
        info = bfgs_o.minimize(x, p);
    }
    // use return_r to get the anisotropy in x
    correct_input23(x, xm, corrections, 1);
    if (p2s_twoIstar)
        tIstar = twoIstar_p2s(expdata->data, M->data, Nchannels);
    else
        tIstar = twoIstar(expdata->data, M->data, Nchannels);
    if (info == 5 || x[0] < 0.) x[0] = -1.;        // for report
    x[1] = xm[1];
    return tIstar;
}


int modelf23(double *param,            // here: [tau gamma r0 rho]
             double *irf,
             double *bg,
             int Nchannels,
             double dt,                // time per channel
             double *corrections,      // [period g l1 l2]
             double *mfunction)        // out: model function in Jordi-girl format
{
#if VERBOSE_FIT2X
    std::cout << "COMPUTE MODEL23" << std::endl;
#endif
    double x[4];
    double tau, gamma, r0, rho,
            period, g, l1, l2,
            taurho, sum_m = 0., tmpf;
    int i, conv_stop;

/************************ Input arguments ***********************/

    tau = param[0];
    gamma = param[1];
    r0 = param[2];
    rho = param[3] * dt / dt;

    period = corrections[0];
    g = corrections[1];
    l1 = corrections[2];
    l2 = corrections[3];
    conv_stop = (int) corrections[4];

#if VERBOSE_FIT2X
    std::cout << "-- tau: " << tau << std::endl;
    std::cout << "-- gamma: " << gamma << std::endl;
    std::cout << "-- r0: " << r0 << std::endl;
    std::cout << "-- rho: " << rho << std::endl;
    std::cout << "-- period: " << period << std::endl;
    std::cout << "-- g-factor: " << g << std::endl;
    std::cout << "-- l1: " << l1 << std::endl;
    std::cout << "-- l2: " << l2 << std::endl;
    std::cout << "-- conv_stop: " << conv_stop << std::endl;
#endif

/************************* Model function ***********************/

    taurho = 1. / (1. / tau + 1. / rho);

    /// vv
    x[0] = 1.;
    x[1] = tau;
    x[2] = r0 * (2. - 3. * l1);
    x[3] = taurho;
    fconv_per_cs(mfunction, x, irf, 2, Nchannels - 1, Nchannels, period, conv_stop, dt);

    /// vh
    x[0] = 1. / g;
    x[2] = 1. / g * r0 * (-1. + 3. * l2);
    fconv_per_cs(mfunction + Nchannels, x, irf + Nchannels, 2, Nchannels - 1, Nchannels, period, conv_stop, dt);

    /// add background
    for (i = 0; i < 2 * Nchannels; i++) sum_m += mfunction[i];
    tmpf = (1. - gamma) / sum_m;
    for (i = 0; i < 2 * Nchannels; i++) mfunction[i] = mfunction[i] * tmpf + bg[i] * gamma;
    return 0;
}

//////////////////////////////////////////// FIT24 ////////////////////////////////////////////

void correct_input24(double *x, double *xm, LVDoubleArray *corrections, int return_r) {
    double r, g, Fp, Fs, l1, l2;

    // correct input parameters (take care of unreasonable values)

    xm[0] = x[0];
    if (xm[0] < 0.001) xm[0] = 0.001;    // tau1 > 0

    xm[2] = x[2];
    if (xm[2] < 0.001) xm[2] = 0.001;    // tau2 > 0

    if (x[3] < 0.) xm[3] = 0.;        // 0 < A2 < 0.999
    else if (x[3] > 0.999) xm[3] = 0.999;
    else xm[3] = x[3];

    if (x[1] < 0.) xm[1] = 0.;        // 0 < gamma < 0.999
    else if (x[1] > 0.999 - xm[3]) xm[1] = 0.999;
    else xm[1] = x[1];

    xm[4] = x[4];
    if (xm[4] < 0.) xm[4] = 0.;        // background > 0

    // anisotropy
    if (return_r) {
        g = corrections->data[1];
        l1 = corrections->data[2];
        l2 = corrections->data[3];
        Fp = (Sp - xm[1] * Bp) / (1. - xm[1]);
        Fs = (Ss - xm[1] * Bs) / (1. - xm[1]);
        r = (Fp - g * Fs) / (Fp * (1. - 3. * l2) + (2. - 3. * l1) * g * Fs);
        x[7] = (Sp - g * Ss) / (Sp * (1. - 3. * l2) + (2. - 3. * l1) * g * Ss);
        x[6] = r;
    }

}


int modelf24(double *param,            // here: [tau1 gamma tau2 A2 offset]
             double *irf,
             double *bg,
             int Nchannels,
             double dt,            // time per channel
             double *corrections,        // [period g l1 l2]
             double *mfunction)        // out: model function in Jordi-girl format

{
    double x[4];
    double tau1, gamma, tau2, A2, offset,
            period,
            sum_m = 0., sum_s = 0.;
    int i, conv_stop;

/************************ Input arguments ***********************/

    tau1 = param[0];
    gamma = param[1];
    tau2 = param[2];
    A2 = param[3];
    offset = param[4] / (double) Nchannels;

    period = corrections[0];
    conv_stop = (int) corrections[4];

/************************* Model function ***********************/

    /// vv
    x[0] = 1. - A2;
    x[1] = tau1;
    x[2] = A2;
    x[3] = tau2;
    fconv_per_cs(mfunction, x, irf, 2, Nchannels - 1, Nchannels, period, conv_stop, dt);

    /// vh
    fconv_per_cs(mfunction + Nchannels, x, irf + Nchannels, 2, Nchannels - 1, Nchannels, period, conv_stop, dt);

    /// add scatter and background

    for (i = 0; i < 2 * Nchannels; i++) {
        sum_m += mfunction[i];
        sum_s += bg[i];
    }
    for (i = 0; i < 2 * Nchannels; i++)
        mfunction[i] = mfunction[i] * (1. - gamma) / sum_m + bg[i] * gamma / sum_s + offset;

    return 0;

}


//////////////////////////////////////////// fit24 ////////////////////////////////////////////

double targetf24(double *x, void *pv) {

    double w, xm[5], Bgamma;
    MParam *p = (MParam *) pv;

    LVI32Array *expdata = *(p->expdata);
    int Nchannels = expdata->length / 2;
    LVDoubleArray *irf = *(p->irf), *bg = *(p->bg),
            *corrections = *(p->corrections), *M = *(p->M);

    correct_input24(x, xm, corrections, 0);
    modelf24(xm, irf->data, bg->data, Nchannels, p->dt, corrections->data, M->data);
    normM_p2s(M->data, Nchannels);

    w = Wcm(expdata->data, M->data, Nchannels);

    if (softbifl & (Bexpected > 0.)) {
        Bgamma = xm[1] * (Sp + Ss);
        w -= Bgamma * log(Bexpected) - loggammaf(Bgamma + 1.);
    }
    return w / Nchannels;

}

double fit24(double *x, short *fixed, MParam *p) {
    // x is:
    // [0] tau1
    // [1] gamma
    // [2] tau2
    // [3] A2 (A1 + A2 = 1)
    // [4] background (offset)
    // [5] BIFL scatter fit? (flag)
    // [6] r Scatter (output only)
    // [7] r Experimental (output only)

    double tIstar, xm[5], B;
    int i, info;

    if (firstcall) init_fact();
    firstcall = 0;
    softbifl = (x[5] < 0.);

    LVI32Array *expdata = *(p->expdata);
    int Nchannels = expdata->length / 2;
    LVDoubleArray *irf = *(p->irf), *bg = *(p->bg),
            *corrections = *(p->corrections), *M = *(p->M);

    // total signal and background

    compute_signal_and_background(p);
    Bexpected = x[1] * (Sp + Ss);

    bfgs bfgs_o(targetf24, 5);

    if (fixed[0]) bfgs_o.fix(0);
    if (fixed[2]) bfgs_o.fix(2);
    if (fixed[3]) bfgs_o.fix(3);
    if (fixed[1] || (x[5] > 0.)) bfgs_o.fix(1);
    if (fixed[4]) bfgs_o.fix(4);

    // fix gamma and offset, try to fit
//  bfgs_o.fix(1);	// gamma
//  bfgs_o.fix(4); 	// background
//  if (!(fixed[0] && fixed[2] && fixed[3])) bfgs_o.minimize(x,p);

    // fit with free gamma and bg
//  if (!fixed[4]) bfgs_o.free(4);
//  if (!fixed[1] && (x[5]<=0.)) bfgs_o.free(1);
//  if ((!fixed[1] && (x[5]<=0.)) || (!fixed[4]))
    info = bfgs_o.minimize(x, p);

    correct_input24(x, xm, corrections, 1);
    modelf24(xm, irf->data, bg->data, Nchannels, p->dt, corrections->data, M->data);
    normM_p2s(M->data, Nchannels);

    tIstar = twoIstar(expdata->data, M->data, Nchannels);

    if (info == 5 || x[0] < 0.) x[0] = -1.;        // for report
    if (info == 5 || x[2] < 0.) x[2] = -1.;
    x[1] = xm[1];
    x[4] = xm[4];

    return tIstar;
}

//////////////////////////////////////////// fit25 ////////////////////////////////////////////

void correct_input25(double* x, double* xm, LVDoubleArray* corrections, int return_r)
{
    double r, g, Fp, Fs, l1, l2;
    // correct input parameters (take care of unreasonable values)
    // here x = [tau gamma r0 rho] + outputs

    xm[0] = x[0];
    penalty = 0.;

    if (x[1]<0.) xm[1] = 0.;		// 0 < gamma < 0.999
    else if (x[1]>0.999) xm[1] = 0.999;
    else xm[1] = x[1];

    xm[2] = x[2];

    // anisotropy
    g = corrections->data[1];
    l1 = corrections->data[2];
    l2 = corrections->data[3];
    Fp = (Sp-xm[1]*Bp)/(1.-xm[1]);
    Fs = (Ss-xm[1]*Bs)/(1.-xm[1]);
    r = (Fp - g*Fs)/(Fp*(1.-3.*l2) + (2.-3.*l1)*g*Fs);

#if VERBOSE_FIT2X
    std::cout << "correct_input25" << std::endl;
    std::cout<< "xm[1]:" << xm[1] << std::endl;
    std::cout<< "Bp:" << Bp << std::endl;
    std::cout<< "Bs:" << Bs << std::endl;
    std::cout<< "Ss:" << Ss << std::endl;
    std::cout<< "Sp:" << Sp << std::endl;
    std::cout<< "Fp:" << Fp << std::endl;
    std::cout<< "Fs:" << Fs << std::endl;
    std::cout<< "l1:" << l1 << std::endl;
    std::cout<< "l2:" << l2 << std::endl;
    std::cout<< "g:" << g << std::endl;
    std::cout<< "gamma:" << x[1] << std::endl;
    std::cout<< "r:" << r << std::endl;
    std::cout<< "rho:" << x[3] << std::endl;
    std::cout<< "tau:" << x[0] << std::endl;
    std::cout<< "r0:" << x[2] << std::endl;
#endif

    if (!fixedrho) {
        xm[3] = x[0]/(x[2]/r-1.);		// rho = tau/(r0/r-1)
        if (xm[3]<1.e-4) xm[3]=1.e-4;	// rho > 0
        x[3] = xm[3];
    }
    else xm[3] = x[3];

    if (return_r) {
        x[8] = (Sp - g*Ss)/(Sp*(1.-3.*l2) + (2.-3.*l1)*g*Ss);
        x[7] = r;
    }

}

double targetf25(double* x, void* pv)
{

    double w, xm[4], Bgamma;
    MParam* p = (MParam*)pv;

    LVI32Array* expdata = *(p->expdata);
    int Nchannels = expdata->length/2;
    LVDoubleArray *irf = *(p->irf), *bg = *(p->bg),
            *corrections = *(p->corrections), *M = *(p->M);

    correct_input25(x, xm, corrections, 0);
    modelf23(xm, irf->data, bg->data, Nchannels, p->dt, corrections->data, M->data);
    normM(M->data, Nchannels);

    if (p2s_twoIstar) w = Wcm_p2s(expdata->data, M->data, Nchannels);
    else w = Wcm(expdata->data, M->data, Nchannels);

    if (softbifl & (Bexpected > 0.)) {
        Bgamma = xm[1]*(Sp+Ss);
        w -= Bgamma*log(Bexpected) - loggammaf(Bgamma+1.);
    }
    return w/Nchannels + penalty;

}


double fit25 (double* x, short* fixed, MParam* p)
{
    // x is:
    // [0] tau1 always fixed
    // [1] tau2 always fixed
    // [2] tau3 always fixed
    // [3] tau4 always fixed
    // [4] gamma
    // [5] r0
    // [6] BIFL scatter fit? (flag)
    // [7] r Scatter (output only)
    // [8] r Experimental (output only)

    double tIstar, tIstarbest = 1.E6, taubest = -1, gammabest = 0., xtmp[9], xm[4], B;
    int i, info;

    if (firstcall) init_fact();
    firstcall = 0;
    softbifl = (x[6]<0.);
    p2s_twoIstar = 1;

    LVI32Array* expdata = *(p->expdata);
    int Nchannels = expdata->length/2;
    LVDoubleArray *irf = *(p->irf), *bg = *(p->bg),
            *corrections = *(p->corrections), *M = *(p->M);

    // total signal and background

    fixedrho = 0;
    Sp = 0.; Ss = 0.; Bp = 0.; Bs = 0.;
    for(i=0; i<Nchannels; i++) {
        Sp += expdata->data[i];
        Bp += bg->data[i];
    }
    for(; i<2*Nchannels; i++) {
        Ss += expdata->data[i];
        Bs += bg->data[i];
    }
    B = Bp+Bs;
    Bp *= (Sp+Ss)/std::max(1., B);
    Bs *= (Sp+Ss)/std::max(1., B);
    Bexpected = x[4]*(Sp+Ss);

    // xtmp: same order as for fit23: [tau gamma r0 rho]
    xtmp[0] = x[0]; xtmp[1] = x[4]; xtmp[2] = x[5]; xtmp[3] = 1.;

    bfgs bfgs_o(targetf25, 4);

    bfgs_o.fix(0);	// tau
    bfgs_o.fix(2); 	// r0
    bfgs_o.fix(3); 	// rho is set in targetf23

    // choose tau
    info = 0;
    for (int i=0; i<4; i++) {
        xtmp[0] = x[i]; xtmp[1] = x[4];
        // fit gamma if unfixed
        if (!fixed[4] && (x[6]<=0.)) bfgs_o.minimize(xtmp,p);

        // calculate 2I*
        correct_input25(xtmp, xm, corrections, 1);
        modelf23(xm, irf->data, bg->data, Nchannels, p->dt, corrections->data, M->data);
        normM(M->data, Nchannels);
        if (p2s_twoIstar) tIstar = twoIstar_p2s(expdata->data, M->data, Nchannels);
        else tIstar = twoIstar(expdata->data, M->data, Nchannels);
//outf << x[i] << '\t' << tIstar << '\t';
#if VERBOSE_FIT2X
        std::cout<< x[i] << "\t" << tIstar << "\t"  << std::endl;
#endif
        if (tIstar < tIstarbest) {
            tIstarbest = tIstar;
            taubest = x[i];
            gammabest = xm[1];
        }
    }
//outf << '\n';

    x[0] = taubest;
    x[4] = gammabest;
    xtmp[0] = x[0]; xtmp[1] = x[4];

    // calculate model function for taubest
    correct_input25(xtmp, xm, corrections, 1);
    modelf23(xm, irf->data, bg->data, Nchannels, p->dt, corrections->data, M->data);
    normM(M->data, Nchannels);
//outf.close();

    x[7] = xtmp[7]; x[8] = xtmp[8];
    return tIstarbest;

}


//////////////////////////////////////////// fit26 ////////////////////////////////////////////

void correct_input26(double* x, double* xm)
{
#if VERBOSE_FIT2X
    std::cout<<"correct_input26"<<std::endl;
#endif
    // correct input parameters (take care of unreasonable values)
    xm[0] = x[0]; // fraction of pattern 1 is between 0 and 1
    if (xm[0]<0.0) {
        xm[0] = 0.0; // tau > 0
        penalty = -x[0];
    }
    else if (xm[0]>1.0) {
        xm[0] = 1.0; // tau > 0
        penalty = x[0]-1.0;
    }
    else penalty = 0.;
#if VERBOSE_FIT2X
    std::cout<<"x[0]: " << x[0] <<std::endl;
    std::cout<<"xm[0]: " << xm[0] <<std::endl;
#endif
}


double targetf26(double* x, void* pv)
{

    double s = 0., xm[1], w, f;
    int i;
    MParam* p = (MParam*)pv;

    LVI32Array* expdata = *(p->expdata);
    int Nchannels = expdata->length;
    LVDoubleArray *irf = *(p->irf), *bg = *(p->bg), *M = *(p->M);

    correct_input26(x, xm);
    f = xm[0];
    // irf is pattern 1, bg is pattern 2
    for(i=0; i<Nchannels; i++)
    {
        M->data[i] = f*irf->data[i] + (1.-f)*bg->data[i];
        s += expdata->data[i];
    }
    for(i=0; i<Nchannels; i++) M->data[i] *= s;

    // divide here Nchannels / 2, because Wcm multiplies Nchannels by two
    w = Wcm(expdata->data, M->data, Nchannels / 2);

    return w/Nchannels + penalty;

}


double fit26 (double* x, short* fixed, MParam* p)
{
    // x is:
    // [0] fraction of pattern 1
    double tIstar, xm[1], f, s = 0., s1 = 0., s2 = 0.;
    int i, info;

    LVI32Array* expdata = *(p->expdata);
    int Nchannels = expdata->length;
    LVDoubleArray *irf = *(p->irf), *bg = *(p->bg), *M = *(p->M);
    for(i=0; i<Nchannels; i++)
    {
        s1 += irf->data[i];
        s2 += bg->data[i];
    }
    s1 = 1./s1; s2 = 1./s2;
    for(i=0; i<Nchannels; i++)
    {
        irf->data[i]*=s1;
        bg->data[i]*=s2;
    }
    bfgs bfgs_o(targetf26, 1);
    info = bfgs_o.minimize(x,p);

    correct_input26(x, xm);
    f = xm[0];

    // irf is pattern 1, bg is pattern 2
    for(i=0; i<Nchannels; i++)
    {
        M->data[i] = f*irf->data[i] + (1.-f)*bg->data[i];
        s += expdata->data[i];
    }
    for(i=0; i<Nchannels; i++) M->data[i] *= s;

    // divide here Nchannels / 2, because twoIstar multiplies Nchannels by two
    tIstar = twoIstar(expdata->data, M->data, Nchannels / 2);
    if (info==5) x[0] = -1.;		// for report
    x[1]=1.-x[0];
    return tIstar;

}
