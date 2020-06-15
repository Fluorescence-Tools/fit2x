// 2I* parameter for Sp + 2Ss

#include <cmath>
#include <iostream>
using namespace std;
#include <time.h>
#include "twoIstar.h"

// factorial
static double fact[150];

// log-likelihood w(C|m) for Cp + 2Cs

double wcm(int C, double mp, double ms)
{
  double s;
  // first term, Cs = 0
  if (C < 150) s = exp(C*log(mp)) / fact[C];
  else s = exp(C*log(mp) - C*log((double)C) + C);
  double w = s;
  double mfactor = ms / (mp * mp);

  int Cp, Csmax = C/2;
  for (int Cs=1; Cs<=Csmax; Cs++) {
    Cp = C - 2*Cs;
    s *= mfactor * (Cp + 2) * (Cp + 1) / (double)Cs;
    w += s;
  }

  return log(w);
}

// overall log-likelihood w(C,M)
double twoIstar_p2s(int* C, double* Mp, double* Ms, int Ndata)
{
  double W = 0., W0 = 0.;
  int nempty = 0;
  for (int i=0; i<Ndata; i++)
    if (C[i]>0) {
      W += wcm(C[i], Mp[i], Ms[i]);
      W0 += wcm(C[i], C[i]*Mp[i]/(Mp[i] + 2.*Ms[i]), C[i]*Ms[i]/(Mp[i] + 2.*Ms[i]));
    }
    else {W += 1.; W0 += 1.;} // nempty++;
  return -2.*(W-W0)/(double)Ndata;
}

int main()
{
  srand((unsigned)time( NULL ));
  init_fact();
  const int Ndata = 256;
  int nt, ntbest, repeats = 100000;
  int* d = new int[Ndata];
  double* Mp = new double[Ndata];
  double* Ms = new double[Ndata];
  double* M = new double[Ndata];
  double* thist = new double[41];
  double s = 0., sm = 0., smean = 0., twoIbest, twoI, taubest;

  for (int i=0; i<41; i++) thist[i] = 0.;

for (int j = 0; j<repeats; j++) {
  for (int i=0; i<Ndata; i++) d[i] = 0.;
  s = 0.;
  double t, dt = 26./Ndata;
  for (int i=0; i<100; i++) {
    t = -4.*log(1.-rand()/(double)RAND_MAX);
    if (t/dt < (double)Ndata) { d[(int)(t/dt)]++; s++; }
  }
  for (int i=0; i<100; i++) {
    t = -4.*log(1.-rand()/(double)RAND_MAX);
    if (t/dt < (double)Ndata) { d[(int)(t/dt)]+=2; s++; }
  }
  smean += s/(double)repeats;
  twoIbest = 100.;
 nt = 0;
for (t=2; t<6; t += 0.1) {
  sm = 0.;
  for (int i=0; i<Ndata; i++) { Mp[i] = exp(-i*dt/t); Ms[i] = exp(-i*dt/t); sm += Mp[i] + Ms[i]; }
  for (int i=0; i<Ndata; i++) { Mp[i] *= s/sm; Ms[i] *= s/sm; M[i] = Mp[i] + 2.*Ms[i];}
  twoI = twoIstar_p2s(d, Mp, Ms, Ndata);
  if (twoI < twoIbest) {twoIbest = twoI; taubest = t; ntbest=nt;}

  nt++;
}
 thist[ntbest]++;
}
  nt = 0;
  for (double t=2; t<6; t += 0.1) {
   cout << thist[nt] << endl;
   nt++;
  }
  // chi2
  double chi2 = 0.;
  for (int i=0; i<Ndata; i++) chi2 += (d[i] - M[i])*(d[i] - M[i])/(Mp[i] + 4.*Ms[i]);
  cout << "chi2 = " << chi2 / (double)Ndata << endl;
  cout << "<s> = " << smean << endl;

  return 1;
}
