//
// Created by Thomas-Otavio Peulen on 6/12/20.
//

#ifndef FIT2X_FIT25_H
#define FIT2X_FIT25_H

/*!
 *
 * @param x
 * @param xm
 * @param corrections
 * @param return_r
 * @return
 */
void correct_input25(double* x, double* xm, LVDoubleArray* corrections, int return_r);

/*!
 *
 * @param x
 * @param pv
 * @return
 */
double targetf25(double* x, void* pv);


/*!
 *
 * @param x
 * @param fixed
 * @param p
 * @return
 */
double fit25 (double* x, short* fixed, MParam* p);

#endif //FIT2X_FIT25_H
