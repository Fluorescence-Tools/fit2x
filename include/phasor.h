//#ifndef FIT2X_PHASOR_H
//#define FIT2X_PHASOR_H
//
//#include <vector>
//#include <cmath>
//#include <algorithm> /* std::max */
//
//#include "tttrlib/TTTR.h" /* TTTR */
//#include "tttrlib/CLSMImage.h" /* CLSMImage */
//
//// namespaces do not work well with SWIG.
//// thus, I use a class with static methods.
//class phasor{
//
//public:
//    /*!
//     * Compute the phasor (g,s) for a selection of micro times
//     *
//     * This function computes the phasor (g,s) for a set of micro times
//     * that are selected out of an vector. The microtimes are selected by a
//     * second vector. The second vector speciefies which indices of the microtime
//     * vector are used to compute the phasor.
//     *
//     * @param micro_times vector of micro times
//     * @param idxs vector of selected indices
//     * @param minimum_number_of_photons
//     * @param frequency the frequency of the phasor
//     * @param g_irf g-value of instrument response phasor
//     * @param s_irf s-value of instrument response phasor
//     * @return vector of length 2: first element g-value, second element s-value
//     */
//    static std::vector<double> compute_phasor(
//            unsigned short *micro_times,
//            std::vector<int> &idxs,
//            double frequency,
//            int minimum_number_of_photons = 1,
//            double g_irf=1.0,
//            double s_irf=0.0
//    );
//
//    /*!
//     * Compute the phasor (g,s) for a all passed micro times
//     *
//     * @param micro_times vector of micro times
//     * @param n_microtimes number of elements in the micro time array
//     * @param frequency the frequency of the phasor
//     * @return vector of length 2: first element g-value, second element s-value
//     */
//    static std::vector<double> compute_phasor_all(
//            unsigned short* microtimes, int n_microtimes,
//            double frequency
//    );
//
//
//    /*!
//     * https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0194578.s001
//     *
//     * @param g_irf g-value of instrument response phasor
//     * @param s_irf s-value of instrument response phasor
//     * @param g_exp
//     * @param s_exp
//     * @return
//     */
//    static double g(
//            double g_irf, double s_irf,
//            double g_exp, double s_exp
//    );
//
//    /*!
//     * https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0194578.s001
//     *
//     * @param g_irf
//     * @param s_irf
//     * @param g_exp
//     * @param s_exp
//     * @return
//     */
//    static double s(
//            double g_irf, double s_irf,
//            double g_exp, double s_exp
//    );
//
//    static void get_phasor_image(
//            float **output, int *dim1, int *dim2, int *dim3, int *dim4,
//            std::shared_ptr<CLSMImage> image,
//            TTTR* tttr_irf = nullptr,
//            double frequency = -1,
//            int minimum_number_of_photons = 30,
//            bool stack_frames = false
//    );
//
//};
//
//
//#endif //FIT2X_PHASOR_H
