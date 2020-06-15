//!  Convolution, scaling, and lamp shift routines
/*!
 * fast and slow convolution routines, with autoscaling  + lamp shift
 * Version under construction!!!
*/
#include <cmath>

#ifndef DELTA_T
#define DELTA_T 0.05
#endif


/*!
 * rescaling -- old version. sum(fit)->sum(decay)
 *
 * @param fit
 * @param decay
 * @param scale
 * @param start
 * @param stop
 */
void rescale(double *fit, double *decay, double *scale, int start, int stop);


/*!
 * rescaling new
 *
 * @param fit
 * @param decay
 * @param w_sq
 * @param scale
 * @param start
 * @param stop
 */
void rescale_w(double *fit, double *decay, double *w_sq, double *scale, int start, int stop);


/*!
 * rescaling -- new version + background. scale = sum(fit*decay/w^2)/sum(fit^2/w^2)
 *
 * @param fit
 * @param decay
 * @param w_sq
 * @param bg
 * @param scale
 * @param start
 * @param stop
 */
void rescale_w_bg(double *fit, double *decay, double *w_sq, double bg, double *scale, int start, int stop);

/*!
 * fast convolution
 *
 * @param fit
 * @param x
 * @param lamp
 * @param numexp
 * @param start
 * @param stop
 */
void fconv(double *fit, double *x, double *lamp, int numexp, int start, int stop);

/*!
 * fast convolution, high repetition rate
 *
 * @param fit
 * @param x
 * @param lamp
 * @param numexp
 * @param start
 * @param stop
 * @param n_points
 * @param period
 */
void fconv_per(
        double *fit, double *x, double *lamp, int numexp, int start, int stop,
        int n_points, double period
);

/*!
 * fast convolution, high repetition rate, with convolution stop
 *
 * @param fit
 * @param x
 * @param lamp
 * @param numexp
 * @param stop
 * @param n_points
 * @param period
 * @param conv_stop
 */
void fconv_per_cs(double *fit, double *x, double *lamp, int numexp, int stop,
                  int n_points, double period, int conv_stop);

/*!
 * fast convolution with reference compound decay
 *
 * @param fit
 * @param x
 * @param lamp
 * @param numexp
 * @param start
 * @param stop
 * @param tauref
 */
void fconv_ref(double *fit, double *x, double *lamp, int numexp, int start, int stop, double tauref);

/*!
 * slow convolution
 *
 * @param fit
 * @param p
 * @param lamp
 * @param start
 * @param stop
 */
void sconv(double *fit, double *p, double *lamp, int start, int stop);

/*!
 * shifting lamp
 *
 * @param lampsh
 * @param lamp
 * @param ts
 * @param n_points
 */
void shift_lamp(double *lampsh, double *lamp, double ts, int n_points);

