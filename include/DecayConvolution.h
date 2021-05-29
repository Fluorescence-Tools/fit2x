#ifndef FIT2X_DECAYCONVOLUTION_H
#define FIT2X_DECAYCONVOLUTION_H

#include <memory>
#include <cmath> /* std::fmod */
#include <iostream>
#include <algorithm> /* std::fill */

#include "fsconv.h"

#include "DecayCurve.h"
#include "DecayModifier.h"
#include "DecayLifetimeHandler.h"


enum convolution_type {
    ConvFastPeriodicTime,
    ConvFastTime,
    ConvFastPeriodic,
    ConvFast,
    ConvFastAVX,
    ConvFastPerAVX
};


class DecayConvolution : public DecayModifier{

public:

    static void compute_corrected_irf(
            DecayCurve *irf,
            DecayCurve *corrected_irf,
            double irf_shift_channels,
            double irf_background_counts
    ) {
#if VERBOSE_FIT2X
        std::clog << "DecayConvolution::compute_corrected_irf:" << std::endl;
        std::clog << "-- irf_shift_channels:" << irf_shift_channels << std::endl;
        std::clog << "-- irf_background_counts:" << irf_background_counts << std::endl;
        std::clog << "-- irf->ptr():" << irf << std::endl;
        std::clog << "-- irf->size():" << irf->size() << std::endl;
        std::clog << "-- corrected_irf->ptr():" << corrected_irf << std::endl;
        std::clog << "-- corrected_irf->size():" << corrected_irf->size() << std::endl;
#endif
        corrected_irf->resize(irf->size());
        for (size_t i = 0; i < irf->size(); i++)
            corrected_irf->_y[i] = std::max(irf->_y[i] - irf_background_counts, 0.0);
        corrected_irf->set_shift(irf_shift_channels);
    }

    //***********************************************//
    //*      STATIC METHODS                         *//
    //***********************************************//
    /*!
     * Compute a mean lifetime by the moments of the decay and the instrument
     * response function.
     *
     * The computed lifetime is the first lifetime determined by the method of
     * moments (Irvin Isenberg, 1973, Biophysical journal).
     *
     * @param irf_histogram
     * @param decay_histogram
     * @param micro_time_resolution
     * @return
     */
    static double compute_mean_lifetime(
            std::vector<double> irf_histogram,
            std::vector<double> decay_histogram,
            double micro_time_resolution
    ){
#if VERBOSE_FIT2X
        std::clog << "DecayConvolution::compute_mean_lifetime" << std::endl;
        std::clog << "-- micro_time_resolution:" << micro_time_resolution << std::endl;
        std::clog << "-- irf_histogram.size():" << irf_histogram.size() << std::endl;
        std::clog << "-- decay_histogram.size():" << decay_histogram.size() << std::endl;
#endif
        double m0_irf = std::accumulate(irf_histogram.begin(), irf_histogram.end(),0.0);
        double m0_data = std::accumulate(decay_histogram.begin(), decay_histogram.end(), 0.0);
#if VERBOSE_FIT2X
        std::clog << "-- m0_irf:" << m0_irf << std::endl;
        std::clog << "-- m0_data:" << m0_data << std::endl;
#endif
        double m1_irf = 0.0;
        double m1_data = 0.0;
        for(size_t i = 0; i<irf_histogram.size(); i++)
            m1_irf += (double) i * irf_histogram[i];
        for(size_t i = 0; i<decay_histogram.size(); i++)
            m1_data += (double) i * decay_histogram[i];
#if VERBOSE_FIT2X
        std::clog << "-- m1_irf:" << m1_irf << std::endl;
        std::clog << "-- m1_data:" << m1_data << std::endl;
#endif
        double g1 = m0_data / m0_irf;
        double g2 = (m1_data - g1 * m1_irf) / m0_irf;
        double tau1 = g2 / g1 * micro_time_resolution;
        return tau1;
    }

private:

    /// Input lifetime spectrum
    DecayLifetimeHandler *lifetime_handler;

    /// Background and shift corrected irf
    DecayCurve *corrected_irf = nullptr;

    /// The time shift of the irf
    double irf_shift_channels = 0.0;

    /// Background counts of irf
    double irf_background_counts = 0.0;

    /// The method used for convolution
    int convolution_method = 0;

    /// The excitation repetition period (usually in nano seconds)
    double excitation_period = std::numeric_limits<double>::max();

    /// Tracks if the computed corrected irf matched the input
    bool corrected_irf_valid = false;

    void convolve_lifetimes(DecayCurve* decay, bool zero_fill = true){
        // Get lifetime spectrum
        auto lh = lifetime_handler->get_lifetime_spectrum();
        auto lt = lh.data(); int ln = lh.size();

        // corrected irf
        auto irfc = get_corrected_irf();

        // convolution range
        auto start = get_start(irfc);
        int stop = get_stop(irfc);
        double dt = decay->get_average_dx();

        // excitation period
        double ex_per = get_excitation_period();

        // convolution method
        int cm = get_convolution_method();

        // Data
        double* my = decay->get_y()->data();
        int nm = (int) decay->get_y()->size();
        double* mx = decay->get_x()->data();

        // IRF
        double* iy = irfc->get_y()->data();
        int ni = (int) irfc->get_y()->size();

        if(nm > 1) {
            if(zero_fill){
                std::fill(my, my + nm, 0.0);
            }
            if (cm == ConvFastPeriodicTime) {
                fconv_per_cs_time_axis(my, nm, mx, nm, iy, ni, lt, ln, start, stop, ex_per);
            } else if (cm == ConvFastTime) {
                fconv_cs_time_axis(my, nm, mx, nm, iy, ni, lt, ln, start, stop);
            } else if (cm == ConvFastPeriodic) {
                fconv_per(my, lt, iy, ln / 2, start, stop, nm, ex_per, dt);
            } else if (cm == ConvFast) {
                fconv(my, lt, iy, ln / 2, start, stop, dt);
            } else if (cm == ConvFastAVX) {
                fconv_avx(my, lt, iy, ln / 2, start, stop, dt);
            } else if (cm == ConvFastPerAVX) {
                fconv_per_avx(my, lt, iy, ln / 2, start, stop, nm, ex_per, dt);
            }
        }
#if VERBOSE_FIT2X
        std::clog << "DecayConvolution::convolve_lifetimes" << std::endl;
        std::clog << "-- convolution_method: " << cm << std::endl;
        std::clog << "-- convolution start, stop: " << start << ", " << stop << std::endl;
        std::clog << "-- histogram spacing (dt): " << dt << std::endl;
        std::clog << "-- excitation_period: " << ex_per << std::endl;
        std::clog << "-- lifetime spectrum: "; for(auto &v: l) std::clog << v << ", "; std::clog << "\n";
        // std::clog << "-- irfc: "; for(auto &v: *(irfc->get_y())) std::clog << v << ", "; std::clog << "\n";
        std::clog << std::endl;
#endif
    }

    void set_data(DecayCurve* v) override {
        DecayModifier::set_data(v);
        auto irf = get_data();
        if(v != nullptr) {
            corrected_irf->set_x(irf->x);
            corrected_irf->set_y(irf->y);
        }
    }

    void update_irf(){
        if(!corrected_irf_valid) {
            compute_corrected_irf(
                    get_irf(), corrected_irf,
                    get_irf_shift_channels(),
                    get_irf_background_counts());
        }
        corrected_irf_valid = true;
    }

public:

    /*================*/
    /* IRF            */
    /*================*/
    void set_irf(DecayCurve *v) {
        set_data(v);
    }

    DecayCurve *get_irf() {
        return get_data();
    }

    /*================*/
    /* Corrected IRF  */
    /*================*/
    void set_irf_shift_channels(double v) {
        double n_channels = std::max(1., (double) get_irf()->size());
        irf_shift_channels = std::fmod(v, n_channels);
        corrected_irf_valid = false;
    }

    double get_irf_shift_channels() const {
        return irf_shift_channels;
    }

    void set_irf_background_counts(double v) {
        irf_background_counts = v;
        corrected_irf_valid = false;
    }

    double get_irf_background_counts() const {
        return irf_background_counts;
    }

    DecayCurve *get_corrected_irf() {
        update_irf();
        return corrected_irf;
    }

    /*================*/
    /* Convolution    */
    /*================*/

    /// The method used for convolution
    /*!
     * 0 - fconv_per_cs_time_axis
     * 1 - fconv_cs_time_axis
     * 2 - fconv_per
     * 3 - fconv
     * 4 - fconv with AVX optimization
     * 5 - fconv_per with AVX optimization
     */
    void set_convolution_method(int v) {
        convolution_method = std::max(0, std::min(5, v));
    }

    int get_convolution_method() const {
        return convolution_method;
    }

    void set_excitation_period(double v) {
        excitation_period = v;
    }

    double get_excitation_period() const {
        return excitation_period;
    }

    double get_mean_lifetime(DecayCurve* decay){
        auto irf = get_corrected_irf()->y;
        auto data = decay->y;
        auto dt = decay->get_average_dx();
        return compute_mean_lifetime(irf, data, dt);
    }

    void set(
            int convolution_method = ConvFastPeriodicTime,
            double excitation_period = 100,
            double irf_shift_channels = 0.0,
            double irf_background_counts = 0
    ) {
        set_irf_background_counts(irf_background_counts);
        set_convolution_method(convolution_method);
        set_excitation_period(excitation_period);
        set_irf_shift_channels(irf_shift_channels);
    }

    DecayConvolution(
            DecayLifetimeHandler *lifetime_handler,
            DecayCurve *instrument_response_function=nullptr,
            int convolution_method = ConvFast,
            double excitation_period = 100,
            double irf_shift_channels = 0.0,
            double irf_background_counts = 0,
            int start = 0,
            int stop = -1,
            bool active = true
    ) : DecayModifier(instrument_response_function, start, stop, active){
#if VERBOSE_FIT2X
        std::clog << "DecayConvolution::DecayConvolution" << std::endl;
#endif
        corrected_irf = new DecayCurve();
        corrected_irf->resize(1);
        default_data->resize(1);
        default_data->y[0] = 1.0;
        this->lifetime_handler = lifetime_handler;
        set(
                convolution_method,
                excitation_period,
                irf_shift_channels,
                irf_background_counts
        );
    }

    ~DecayConvolution() override {
        delete corrected_irf;
    }

    void add(DecayCurve* out){
        if(is_active()){
            // resize output to IRF
            out->resize(get_data()->size(), 0.0);
            convolve_lifetimes(out);
        }
    }

};

#endif //FIT2X_DECAYCONVOLUTION_H
