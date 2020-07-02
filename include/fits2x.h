#define FIT2X_VERSION                "0.0.4"
//#define VERBOSE 0

#include <iostream>
#include <cmath>
#include <algorithm>

#include "i_lbfgs.h"
#include "fsconv.h"
#include "lvarrays.h"
#include "twoIstar.h"


/*!
 * Computes the total number of photons in the parallel and perpendicular
 * detection channel for the background and the measured signal. The
 * computed number of photons are stored in the static variables
 * Sp, Ss, Bp, Bs.
 *
 * @param p[in] a pointer to a MParam object
 */
void compute_signal_and_background(MParam* p);


/*!
 * Normalizes the number of photons in the entire model function to the
 * number of experimental photons.
 *
 * Here, the Number of experimental photons is Sp + Ss (signal in parallel and
 * perpendicular). Sp and Ss are global variables that can be computed by
 * `compute_signal_and_background`.
 *
 * @param M[in,out] array containing the model function in Jordi format
 * @param Nchannels[in] number of channels in the experiment(half length of
 * M array)
 */
void normM(double* M, int Nchannels);


/*!
 * Normalizes a model function (that is already normalized to a unit area) to
 * the total number of photons in parallel and perpendicular,
 *
 * @param M[in,out] array containing the model function in Jordi format
 * @param s[in] a scaling factor by which the model function is divided.
 * @param Nchannels[in] the number of channels in the model function (half length of
 * M array)
 */
void normM(double* M, double s, int Nchannels);


/*!
 * Normalizes the number of photons in the model function for Ss and Sp
 * individually to the number of experimental photons in Ss and Sp.
 *
 * Here, the number of experimental photons are global variables that can be
 * computed by `compute_signal_and_background`.
 *
 * @param M array[in,out] containing the model function in Jordi format
 * @param Nchannels[in] number of channels in the experiment (half length of
 * M array)
 */
void normM_p2s(double* M, int Nchannels);

#include "fit23.h"
#include "fit24.h"
#include "fit25.h"
#include "fit26.h"

