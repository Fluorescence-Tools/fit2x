//
// Created by Thomas-Otavio Peulen on 6/12/20.
//

#ifndef FIT2X_FIT24_H
#define FIT2X_FIT24_H

/*!
 * Bi-exponential model function
 * @param param
 * @param irf
 * @param bg
 * @param Nchannels
 * @param dt
 * @param corrections
 * @param mfunction
 * @return
 */
int modelf24(double* param, 			// here: [tau1 gamma tau2 A2 offset]
             double* irf,
             double* bg,
             int Nchannels,
             double dt,			// time per channel
             double* corrections,		// [period g l1 l2]
             double* mfunction);		// out: model function in Jordi-girl format


/*!
*
* @param x
* @param pv
* @return
*/
double targetf24(double *x, void *pv);


/*!
 *
 * @param x
 * @param fixed
 * @param p
 * @return
 */
double fit24(double *x, short *fixed, MParam *p);


/*!
 *
 * @param x
 * @param xm
 * @param corrections
 * @param return_r
 * @return
 */
void correct_input24(double *x, double *xm, LVDoubleArray *corrections, int return_r);


#endif //FIT2X_FIT24_H
