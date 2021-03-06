#ifndef FIT2X_FIT26_H
#define FIT2X_FIT26_H

/*!
 * Correct input for fit 26
 *
 * Constrains the fraction x1 of the first pattern to (0 < x1 < 1).
 *
 * @param x[in] x[0] fraction of first pattern
 * @param xm[out] xm[0] corrected fraction of first pattern
 */
void correct_input26(double* x, double* xm);

double targetf26(double* x, void* pv);

/*!
 * Pattern-fit
 *
 * Fits the fraction of a mixture of two patterns
 *
 * The two patterns are set by the attributes irf and bg of the MParam
 * structure.
 *
 * @param x [0] fraction of pattern 1
 * @param fixed not used
 * @param p an instance of MParam that contains the patterns. The fist pattern is
 * contained in the instrument response function array, the second in the background,
 * array, the experimental data is in the array expdata.
 * @return
 */
double fit26(double* x, short* fixed, MParam* p);

#endif //FIT2X_FIT26_H
