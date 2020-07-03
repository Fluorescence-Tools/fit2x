#ifndef TTTRLIB_STATISTICS_H
#define TTTRLIB_STATISTICS_H


#include <cmath>
#include <vector>
#include <numeric>
#include <string>
#include <algorithm>
#include <iostream>



/*!
 * Initialize an array containing pre-computed logratithms
 */
void init_fact();

/*!
 * Approximation of log(gamma function). See wikipedia
 *
 * https://en.wikipedia.org/wiki/Gamma_function#The_log-gamma_function
 *
 * @param t input of the gamma function
 * @return approximation of the logarithm of the gamma function
 */
double loggammaf(double t);

/*!
 * log-likelihood w(C|m) for Cp + 2Cs
 *
 * @param C number of counts in channel
 * @param m model function
 * @return log-likelihood w(C|m) for Cp + 2Cs
 */
double wcm(int C, double m);

/*!
 * Compute the -log-likelihood for Cp + 2Cs of a single micro time channel.
 *
 * Compute score of model counts in a parallel and perpendicular detection channel
 * and the experimental counts for a micro time channel.
 *
 * This function computes a score for the experimental counts (C) in a channel
 * where the experimental counts were computed by the sum of the counts in the
 * parallel (P) and the perpendicular (S) channel by the equation C = P + 2 S.
 *
 * This function considers that the number of counts C = P + 2S is not Poissonian.
 * The score relates to a maximum likelihood function.
 *
 * @param C number of experimental counts (P + 2 S) in a micro time channel
 * @param mp number of counts of the model in parallel detection channel
 * @param ms number of counts of the model in the perpendicular detection channel
 * @return
 */
double wcm_p2s(int C, double mp, double ms);

/*!
 * Compute the overall -log-likelihood for Cp + 2Cs for all micro time channels
 *
 * @param C array of experimental counts in Jordi format
 * @param M array model function in Jordi format
 * @param Nchannels number of micro time channels in parallel and perpendicular
 * (half the number of elements in C and M).
 * @return -log-likelihood for Cp + 2Cs for all micro time channels
 */
double Wcm_p2s(int* C, double* M, int Nchannels);

/*!
 * Compute overall 2I* for Cp + 2Cs
 *
 * This function computes the overall 2I* for the model function Cp + 2Cs that
 * is computed by parallel signal (Cp) and the perpendicular signal (Cs). For
 * the definition of function 2I* see "An Experimental Comparison of the
 * Maximum Likelihood Estimation and Nonlinear Least-Squares Fluorescence Lifetime
 * Analysis of Single Molecules, Michael Maus, Mircea Cotlet, Johan Hofkens,
 * Thomas Gensch, Frans C. De Schryver, J. Schaffer, and C. A. M. Seidel, Anal.
 * Chem. 2001, 73, 9, 2078–2086".
 *
 * @param C array of experimental counts in Jordi format
 * @param M array model function in Jordi format
 * @param Nchannels number of micro time channels in parallel and perpendicular
 * (half the number of elements in C and M).
 * @return 2I* for Cp + 2Cs
 */
double twoIstar_p2s(int* C, double* M, int Nchannels);

/*!
 * Compute overall 2I* for Cp & Cs
 *
 * This function computes 2I* for Cp and Cs. Cp and Cs are the model signals in
 * the parallel and perpendicular channel. Contrary to twoIstar_p2s the overall
 * 2I* is the sum of 2I* for Cp and Cs.
 *
 * @param C array of experimental counts in Jordi format
 * @param M array model function in Jordi format
 * @param Nchannels number of micro time channels in parallel and perpendicular
 * (half the number of elements in C and M).
 * @return 2I* for Cp & Cs
 */
double twoIstar(int* C, double* M, int Nchannels);

/*!
 * Compute overall -log-likelihood for Cp & Cs
 *
 * @param C array of experimental counts in Jordi format
 * @param M array model function in Jordi format
 * @param Nchannels number of micro time channels in parallel and perpendicular
 * (half the number of elements in C and M).
 * @return -log-likelihood for Cp & Cs
 */
double Wcm(int* C, double* M, int Nchannels);

namespace statistics{

    /*!
     * Different chi2 measures for counting data:
     *
     * https://arxiv.org/pdf/1903.07185.pdf
     *
     * @param data
     * @param model
     * @param x_min
     * @param x_max
     * @param type
     * @return
     */
    double chi2_counting(
            std::vector<double> &data,
            std::vector<double> &model,
            int x_min = -1,
            int x_max = -1,
            std::string type="neyman"
    );
}

#endif //TTTRLIB_STATISTICS_H
