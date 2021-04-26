#ifndef FIT2X_DECAYCONVOLUTION_H
#define FIT2X_DECAYCONVOLUTION_H

#include <memory>
#include <cmath> /* std::fmod */
#include <iostream>

#include "fsconv.h"

#include "DecayCurve.h"
#include "DecayLifetimeSpectrum.h"


enum convolution_type{
    ConvFastPeriodicTime,
    ConvFastTime,
    ConvFastPeriodic,
    ConvFast,
    ConvFastAVX,
    ConvFastPerAVX
};


class DecayConvolution{

    friend class DecayCurve;
    friend class Decay;

public:

    static void compute_corrected_irf(
            DecayCurve& irf,
            DecayCurve& corrected_irf,
            double irf_shift_channels,
            double irf_background_counts
    ){
        corrected_irf.resize(irf.size());
        for (int i = 0; i < irf.size(); i++)
            corrected_irf.y[i] = std::max(irf.y[i] - irf_background_counts, 0.0);
        corrected_irf.set_shift(irf_shift_channels);
    }

private:

    DecayLifetimeSpectrum* lifetime_spectrum;

    /// Output Decay curve
    DecayCurve decay;

    /// Tracks if the decay is valid
    bool decay_valid = false;

    /// Input instrument response function
    DecayCurve* irf;

    /// Background and shift corrected irf
    DecayCurve corrected_irf;

    /// Tracks if the corrected irf is valid
    bool corrected_irf_valid = false;

    /// If set to true uses the background corrected IRF as scatter
    bool use_corrected_irf_as_scatter = true;

    /// The time shift of the irf
    double irf_shift_channels = 0.0;

    /// Background counts of irf
    double irf_background_counts = 0.0;

    /// The fraction of the irf (scattered light) in the model decay
    double scatter_fraction = 0.0;

    /// Convolution range
    int convolution_start = 0;
    int convolution_stop = -1;

    /// The method used for convolution
    int convolution_method = 0;

    /// The excitation repetition period (usually in nano seconds)
    double excitation_period = std::numeric_limits<double>::max();

    void update_corrected_irf(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::update_corrected_irf" << std::endl;
#endif
        compute_corrected_irf(
                *irf, corrected_irf,
                get_irf_shift_channels(),
                get_irf_background_counts());
        corrected_irf_valid = true;
    }

    void add_scatter();

    void convolve_lifetimes();

public:

    void update_decay(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::update_decay" << std::endl;
#endif
        convolve_lifetimes();
        add_scatter();
    }

    /// The method used for convolution
    /*!
     * 0 - fconv_per_cs_time_axis
     * 1 - fconv_cs_time_axis
     * 2 - fconv_per
     * 3 - fconv
     * 4 - fconv with AVX optimization
     * 5 - fconv_per with AVX optimization
     */
    void set_convolution_method(int v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_convolution_method" << std::endl;
#endif
        convolution_method = v;
        decay_valid = false;
    }

    int get_convolution_method(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_convolution_method" << std::endl;
#endif
        return convolution_method;
    }

    void set_excitation_period(double v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_excitation_period" << std::endl;
#endif
        excitation_period = v;
        decay_valid = false;
    }

    double get_excitation_period(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_excitation_period" << std::endl;
#endif
        return excitation_period;
    }

    void set_irf_shift_channels(double v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_irf_shift_channels" << std::endl;
#endif
        corrected_irf_valid = false;
        decay_valid = false;
        double n_channels = std::max(1., (double) irf->size());
        irf_shift_channels = std::fmod(v, n_channels);
    }

    double get_irf_shift_channels(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_irf_shift_channels" << std::endl;
#endif
        return irf_shift_channels;
    }

    void set_scatter_fraction(double v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_scatter_fraction" << std::endl;
#endif
        v = std::min(1., std::max(0.0, v));
        scatter_fraction = v;
        decay_valid = false;
    }

    double get_scatter_fraction(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_scatter_fraction" << std::endl;
#endif
        return scatter_fraction;
    }

    void set_convolution_start(int v) {
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_convolution_start" << std::endl;
#endif
        convolution_start = v;
        decay_valid = false;
    }

    int get_convolution_start() const {
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_convolution_start" << std::endl;
#endif
        int start = convolution_start;
        int nmax = irf->size();
        if(start < 0){
            start = nmax - start;
        }
        return std::min(start, nmax);
    }

    void set_convolution_stop(int v) {
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_convolution_stop" << std::endl;
#endif
        convolution_stop = v;
        decay_valid = false;
    }

    int get_convolution_stop() const {
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_convolution_stop" << std::endl;
#endif
        int stop = convolution_stop;
        int nmax = irf->size();
        if(stop < 0){
            stop = std::min(nmax, std::max(0, nmax + stop));
        } else if (stop == 0){
            stop = nmax;
        }
        return std::min(stop, nmax);
    }

    void set_irf_background_counts(double v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_irf_background_counts" << std::endl;
#endif
        corrected_irf_valid = false;
        decay_valid = false;
        irf_background_counts = v;
    }

    double get_irf_background_counts(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_irf_background_counts" << std::endl;
#endif
        return irf_background_counts;
    }

    void set_irf(DecayCurve* v) {
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_irf" << std::endl;
#endif
        irf = v;
        decay.resize(v->size());
        corrected_irf.resize(v->size());
        decay.set_x(irf->x);
        corrected_irf.set_x(irf->x);
        corrected_irf.set_y(irf->y);
        corrected_irf_valid = false;
        decay_valid = false;
    }

    void set_irf(double *input, int n_input){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_irf" << std::endl;
#endif
        decay.resize(n_input);
        irf->set_y(input, n_input);
        decay.resize(n_input);
        corrected_irf.resize(n_input);
        corrected_irf.set_y(irf->y);
        corrected_irf_valid = false;
        decay_valid = false;
    }

    DecayCurve* get_irf(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_irf" << std::endl;
#endif
        return irf;
    }

    void get_irf(double **output_view, int *n_output){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_irf" << std::endl;
#endif
        irf->get_y(output_view, n_output);
    }

    DecayCurve* get_corrected_irf(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_corrected_irf" << std::endl;
#endif
        if(!corrected_irf_valid){
            update_corrected_irf();
        }
        return &corrected_irf;
    }

    void get_corrected_irf(double **output_view, int *n_output){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_corrected_irf" << std::endl;
#endif
        corrected_irf.get_y(output_view, n_output);
    }

    void set_lifetime_spectrum(DecayLifetimeSpectrum* v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_lifetime_spectrum" << std::endl;
#endif
        lifetime_spectrum = v;
        decay_valid = false;
    }

    DecayLifetimeSpectrum* get_lifetime_spectrum(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_lifetime_spectrum" << std::endl;
#endif
        return lifetime_spectrum;
    }

    void set_use_corrected_irf_as_scatter(bool v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set_use_corrected_irf_as_scatter" << std::endl;
#endif
        use_corrected_irf_as_scatter = v;
        decay_valid = false;
    }

    bool get_use_corrected_irf_as_scatter(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_use_corrected_irf_as_scatter" << std::endl;
#endif
        return use_corrected_irf_as_scatter;
    }

    DecayCurve* get_decay(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::get_decay" << std::endl;
#endif
        if(!decay_valid){
            update_decay();
        }
        return &decay;
    }

    void set(
            DecayCurve* instrument_response_function,
            DecayLifetimeSpectrum* lifetime_spectrum,
            std::vector<int> convolution_range = std::vector<int>({0, -1}),
            bool use_corrected_irf_as_scatter = true,
            double scatter_fraction = 0.0,
            int convolution_method = ConvFastPeriodicTime,
            double excitation_period = 100,
            double irf_shift_channels = 0.0,
            double irf_background_counts = 0
    ){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::set" << std::endl;
#endif
        set_irf(instrument_response_function);
        set_lifetime_spectrum(lifetime_spectrum);
        set_irf_background_counts(irf_background_counts);
        set_convolution_start(convolution_range[0]);
        set_convolution_stop(convolution_range[1]);
        set_use_corrected_irf_as_scatter(use_corrected_irf_as_scatter);
        set_scatter_fraction(scatter_fraction);
        set_convolution_method(convolution_method);
        set_excitation_period(excitation_period);
        set_irf_shift_channels(irf_shift_channels);
    }

    DecayConvolution(
            DecayCurve* instrument_response_function,
            DecayLifetimeSpectrum* lifetime_spectrum,
            std::vector<int> convolution_range = std::vector<int>({0, -1}),
            bool use_corrected_irf_as_scatter = true,
            double scatter_fraction = 0.0,
            int convolution_method = ConvFastPeriodicTime,
            double excitation_period = 100,
            double irf_shift_channels = 0.0,
            double irf_background_counts = 0
            ){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayConvolution::DecayConvolution" << std::endl;
#endif
        set(
                instrument_response_function,
                lifetime_spectrum,
                convolution_range,
                use_corrected_irf_as_scatter,
                scatter_fraction,
                convolution_method,
                excitation_period,
                irf_shift_channels,
                irf_background_counts
        );
    }

};

#endif //FIT2X_DECAYCONVOLUTION_H
