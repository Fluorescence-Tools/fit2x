//
// Created by Thomas-Otavio Peulen on 6/12/20.
//

#ifndef FIT2X_FIT26_H
#define FIT2X_FIT26_H

void correct_input26(double* x, double* xm);

double targetf26(double* x, void* pv);

double fit26 (double* x, short* fixed, MParam* p);

#endif //FIT2X_FIT26_H
