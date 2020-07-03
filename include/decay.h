#ifndef TTTRLIB_DECAY_H
#define TTTRLIB_DECAY_H

#include <cmath>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <numeric>
#include <algorithm>
#include <memory> /* shared_ptr */

#include "omp.h"
#include "tttrlib/tttr.h"
#include "statistics.h"
#include "fsconv.h"



class Decay {

private:

    /// Stores the weighted residuals
    std::vector<double> _weighted_residuals = std::vector<double>();

    /// Proportional to the area of the model. If negative, the area is scaled
    /// to the data.
    double _total_area = -1;

    /// The constant offset in the model
    double _constant_background = 0.0;

    /// The fraction of the irf (scattered light) in the model decay
    double _areal_scatter_fraction = 0.0;

    /// The time shift of the irf
    double _irf_shift_channels = 0.0;

    /// Background counts of irf. This number is sutracted from the counts
    /// in each channel before convolution.
    double _irf_background_counts = 0.0;

    /// The repetiotion period (usually in nano seconds)
    double _period = 100.;

    /// The instrument response function
    std::vector<double> _irf = std::vector<double>();

    /// The experimental data (used for scaling)
    std::vector<double> _data = std::vector<double>();

    /// The weights of the experimental data (used for weighted residuals)
    std::vector<double> _weights = std::vector<double>();

    /// The squared weights of the experimental data (used for scaling)
    std::vector<double> _sq_weights = std::vector<double>();

    /// The time axis of the data (used in convolution)
    std::vector<double> _time_axis = std::vector<double>();

    /// The used to store the model function
    std::vector<double> _model_function = std::vector<double>();

    /// The lifetime spectrum
    std::vector<double> _lifetime_spectrum = std::vector<double>();

    /// Set to valid if the output decay matched the input
    bool _is_valid = false;

    /// If set to true will add pile up to the model function
    bool _correct_pile_up = false;

    /// Beginning index of the comvolution
    int _convolution_start = 0;

    /// Stop index of the convolution
    int _convolution_stop = -1;

    /// Amplitude threshold used to discriminate low populated lifetimes
    double _amplitude_threshold = 1e10;

    /// If set to true discriminates low populated lifetimes in convolution.
    bool _use_amplitude_threshold = false;


public:

    /// The range on which the scoring function is evaluated
    int x_min = -1;
    int x_max = -1;

    /// If set to true take abs values of lifetime spectrum
    bool abs_lifetime_spectrum = true;

    /*!
     * Computes a fluorescence decay for a lifetime spectrum.
     *
     * The lifetime spectrum is a interleaved array of amplitudes and fluorescence
     * lifetimes. The fluorescence decay for the lifetime spectrum is computed.
     * The computed decay is convolved with a instrument response function (IRF).
     * Before convolution the IRF is corrected for a constant offset. The convolution
     * considers periodic excitation. The IRF is shifted by a specified value of
     * micro time channels. After convolution a constant fraction of scattered
     * light is added to the computed decay and the decay is scaled to the number
     * of photons in an experimental fluorescence decay. Finally, a constant
     * background is added. Optionally, pile-up is added to the computed fluorescence
     * decay. Before the fluorescence decay is computed an amplitude threshold can
     * be applied to the fluorescence lifetime spectrum to discriminate fluorescence
     * lifetimes with small amplitudes.
     *
     * @param model_function[out] The output array that will contain the computed
     * fluorescence decay
     * @param n_model_function[out] The number of points in the model function
     * @param data[in] The data for which the model function is computed. The data
     * will be used to scaled the computed decay.
     * @param n_data[in] The number of data points
     * @param squared_weights[in] The squared weights of the data points. The data
     * weights are used to scale the model function to the data (usually the data weights is Poissonian)
     * @param n_weights[in] The number of weights
     * @param time_axis[in] The time axiss for which the fluorescence decay is computed
     * for
     * @param n_time_axis[in] The number of points in the time axis
     * @param instrument_response_function[in] The instrument response function
     * used for convolution
     * @param n_instrument_response_function[in] The number of points in the
     * instrument response function
     * @param lifetime_spectrum[in] The lifetime spectrum. A lifetime specturm is an
     * interleaved array of amplitudes and fluorescence lifetimes, .e.g,
     * (amplitude1, lifetime1, amplitude2, lifetime2, ...)
     * @param n_lifetime_spectrum[in] The number of points in the fluorescence lifetime
     * spectrum
     * @param start[in] The start of the convolution
     * @param stop[in] The stop of the convoution. The decay will be computed in the
     * range [start, stop]
     * @param irf_background_counts[in] The background counts of the instrument
     * response function. This number will be subtracted from the IRF before
     * convolution.
     * @param irf_shift_channels[in] The number of micro time channels the IRF will be
     * shifted before the fluorescence lifetimes are convoluted with the IRF.
     * @param scatter_areal_fraction[in] The fraction (integrated fraction), i.e.,
     * the area the scattered light will have in the computed decay.
     * @param excitation_period[in] The excitation period (in units of the fluorescence
     * lifetime, usually nanoseconds) that was used to excite the sample
     * @param constant_background A constant offset that is added to the fluorescence
     * decay.
     * @param total_area[in] the area to which the model fluorescence decay is scaled
     * to. If this value is negative (default value is -1) the model fluorescence
     * decay is scaled to the data in the range defined by the parameters start,
     * stop
     * @param use_amplitude_threshold[in] if set to true (default is false) a
     * threshold will be applied. If this parameter is set to true, fluorescence
     * lifetimes with a amplitude that is smaller than a threshold will be not
     * considered.
     * @param amplitude_threshold The threshold that is used to discriminate
     * fluorescence lifetimes that are smaller.
     * @param pile_up if set to true (default is false) pile up will be added
     * to the model function.
     * @param instrument_dead_time the dead time of the instrument (used for
     * pile up (in units of the lifetime, usually nano seconds)
     * @param acquisition_time the total time the acquisition of the decay in
     * seconds.
     * @param add_corrected_irf if set to true (default is false) the background
     * corrected irf will be added as scatter to the decay. If this is false the
     * irf prior to a background and shift corrected irf is added as a scatter
     * fraction.
     * @param scale_model_to_area if set to true (default is true) the model is
     * either scaled to the area s provided by total_area (if total_area is larger
     * then zero) or to the total number of experimental counts.
     */
    static void compute_decay(
            double *model_function, int n_model_function,
            double *data, int n_data,
            double *squared_weights, int n_weights,
            double *time_axis, int n_time_axis,
            double *instrument_response_function, int n_instrument_response_function,
            double *lifetime_spectrum, int n_lifetime_spectrum,
            int start = 0, int stop = -1,
            double irf_background_counts = 0.0,
            double irf_shift_channels = 0.0,
            double scatter_areal_fraction = 0.0,
            double excitation_period = 1000.0,
            double constant_background = 0.0,
            double total_area = -1,
            bool use_amplitude_threshold = false,
            double amplitude_threshold = 1e10,
            bool pile_up = false,
            double instrument_dead_time=120.0,
            double acquisition_time=100000.0,
            bool add_corrected_irf = false,
            bool scale_model_to_area = true
    );

protected:

    bool is_valid() const {
        return _is_valid;
    }

public:

    void set_data(double *input, int n_input) {
        _is_valid = false;
        _data.assign(input, input + n_input);
    }

    void get_data(double **output_view, int *n_output) {
        *output_view = _data.data();
        *n_output = _data.size();
    }

    void set_use_amplitude_threshold(bool v) {
        _is_valid = false;
        _use_amplitude_threshold = v;
    }

    bool get_use_amplitude_threshold() const {
        return _use_amplitude_threshold;
    }

    void set_amplitude_threshold(double v) {
        _is_valid = false;
        _amplitude_threshold = v;
    }

    double get_amplitude_threshold() const {
        return _amplitude_threshold;
    }

    void set_total_area(double v) {
        _is_valid = false;
        _total_area = v;
    }

    double get_total_area() const {
        return _total_area;
    }

    void set_period(double v) {
        _is_valid = false;
        _period = v;
    }

    void set_score_range(int vmin, int vmax){
        x_min = vmin; x_max = vmax;
    }

    double get_period() const {
        return _period;
    }

    void set_irf_shift_channels(double v) {
        _is_valid = false;
        _irf_shift_channels = std::fmod(v, (double)_irf.size());
    }

    double get_irf_shift_channels() const {
        return _irf_shift_channels;
    }

    void set_areal_scatter_fraction(double v) {
        _is_valid = false;
        _areal_scatter_fraction = std::min(1., std::max(0.0, v));
    }

    double get_areal_scatter_fraction() const {
        return _areal_scatter_fraction;
    }

    void set_constant_background(double v) {
        _is_valid = false;
        _constant_background = v;
    }

    double get_constant_background() const {
        return _constant_background;
    }

    void set_convolution_start(int v) {
        _is_valid = false;
        _convolution_start = v;
    }

    int get_convolution_start() const {
        return _convolution_start;
    }

    void set_convolution_stop(int v) {
        _is_valid = false;
        _convolution_stop = v;
    }

    int get_convolution_stop() const {
        return _convolution_stop;
    }

    void set_correct_pile_up(bool v) {
        _is_valid = false;
        _correct_pile_up = v;
    }

    bool get_correct_pile_up() const {
        return _correct_pile_up;
    }

    void set_irf(double *input, int n_input) {
        _is_valid = false;
        _irf.clear();
        _irf.assign(input, input + n_input);
    }

    void get_irf(double **output_view, int *n_output) {
        *output_view = _irf.data();
        *n_output = _irf.size();
    }

    void get_model(double **output_view, int *n_output) {
        if (!_is_valid) {
            evaluate();
        }
        *output_view = _model_function.data();
        *n_output = _model_function.size();
    }

    void set_lifetime_spectrum(double *input, int n_input = -1) {
        _is_valid = false;
        if(n_input < 0) n_input = _lifetime_spectrum.size();
        _lifetime_spectrum.resize(n_input);
        for (int i = 0; i < n_input; i++) {
            _lifetime_spectrum[i] = input[i];
        }
    }

    void get_lifetime_spectrum(double **output_view, int *n_output) {
        *output_view = _lifetime_spectrum.data();
        *n_output = _lifetime_spectrum.size();
    }

    void set_weights(double *input, int n_input) {
        _is_valid = false;
        _weights.clear();
        _weights.assign(input, input + n_input);
    }

    void get_weights(double **output_view, int *n_output) {
        *output_view = _weights.data();
        *n_output = _weights.size();
    }

    void set_time_axis(double *input, int n_input) {
        _is_valid = false;
        _time_axis.assign(input, input + n_input);
    }

    void get_time_axis(double **output_view, int *n_output) {
        *output_view = _time_axis.data();
        *n_output = _time_axis.size();
    }

    void set_irf_background_counts(
            double irf_background_counts
    ) {
        _irf_background_counts = irf_background_counts;
        _is_valid = false;
    }

    double get_irf_background_counts() const {
        return _irf_background_counts;
    }

    /*!
     *
     * @param tttr_data pointer to TTTR object that is used to construct a decay
     * histogram
     * @param micro_time_coarsening an (optional) integer by which the micro times
     * are divided to coarsen the time axis (default is 1)
     * @param decay_histogram the data to which the decay is fitted
     * @param time_axis the time axis that belongs to the data
     * @param dt the spacing between the points in the time axis. This optional
     * parameter is used to compute a time axis if not time axis was provided by
     * the parameter time_axis
     * @param weights the weights of the data points. If the weights are not provided
     * (nullptr / None) the weights are computed assuming Poisson noise.
     * @param instrument_response_function The instrument response function (IRF)
     * that is used for convolution. If no IRF is provided
     * @param start The start index in the IRF used for convolution. Points in the
     * IRF before the start index are not used for convolution.
     * @param stop The stop index in the IRF used for convolution. Points beyond the
     * stop index are not convolved.
     * @param use_amplitude_threshold If this is set to true (default value is true)
     * the values that are smaller then a specified threshold are omitted
     * @param amplitude_threshold The amplitude threshold that is used if the
     * parameter use_amplitude_threshold is set to true (the default value is 1e10)
     * @param correct_pile_up If this is set to true (the default value is false)
     * the convolved model function is 'piled up' to match pile up artifacts in the
     * data.
     * @param excitation_period the repetition period, .i.e, the time between subsequent
     * excitation pulses.
     */
    Decay(
            std::shared_ptr<TTTR> tttr_data = nullptr,
            std::shared_ptr<TTTR> tttr_irf = nullptr,
            int micro_time_coarsening = 1,
            std::vector<double> decay_histogram = std::vector<double>(),
            std::vector<double> time_axis = std::vector<double>(),
            std::vector<double> weights = std::vector<double>(),
            std::vector<double> instrument_response_function = std::vector<double>(),
            int start = 0, int stop = -1,
            bool use_amplitude_threshold = false,
            double amplitude_threshold = 1e10,
            bool correct_pile_up = false,
            double excitation_period = 100.0,
            double dt = 1.0,
            double irf_background_counts = 0.0,
            double areal_scatter_fraction = 0.0,
            double constant_background = 0.0,
            int x_min = -1,
            int x_max = -1,
            std::vector<double> lifetime_spectrum = std::vector<double>()
    ) {
        _convolution_start = std::max(0, start);
        _convolution_stop = std::min(stop, (int) decay_histogram.size());
        _amplitude_threshold = amplitude_threshold;
        _use_amplitude_threshold = use_amplitude_threshold;
        _correct_pile_up = correct_pile_up;
        _period = (tttr_data == nullptr) ?
                  excitation_period :
                  tttr_data->get_header().macro_time_resolution;
        _irf_background_counts = irf_background_counts;
        _areal_scatter_fraction = areal_scatter_fraction;
        _constant_background = constant_background;
        this->x_min = x_min;
        this->x_max = x_max;
        set_lifetime_spectrum(lifetime_spectrum.data(), lifetime_spectrum.size());
        // set data
        if (tttr_data != nullptr) {
            double *hist; int n_hist;
            double *time; int n_time;
            TTTR::compute_microtime_histogram(
                    tttr_data.get(),
                    &hist, &n_hist,
                    &time, &n_time,
                    micro_time_coarsening
            );
            set_data(hist, n_hist);
            set_time_axis(time, n_time);
            free(hist); free(time);
        } else {
            set_data(decay_histogram.data(), decay_histogram.size());
        }
        // set time axis
        if(time_axis.empty()) {
#if VERBOSE
            std::clog << "-- Filling time axis, dt: " << dt << std::endl;
#endif
            _time_axis.clear();
            _time_axis.reserve(_data.size());
            for (int i = 0; i < _data.size(); i++) _time_axis.emplace_back(i * dt);
        } else{
            set_time_axis(time_axis.data(), time_axis.size());
        }
        // set irf
        if (tttr_irf != nullptr) {
#if VERBOSE
            std::clog << "-- Setting IRF from TTTR..." << std::endl;
#endif
            double *hist; int n_hist;
            double *time; int n_time;
            TTTR::compute_microtime_histogram(
                    tttr_irf.get(),
                    &hist, &n_hist,
                    &time, &n_time,
                    micro_time_coarsening
            );
            set_irf(hist, n_hist);
            free(time); free(hist);
        } else {
#if VERBOSE
            std::clog << "-- Setting IRF from argument" << std::endl;
#endif
            if (instrument_response_function.empty() && !_data.empty()) {
#if VERBOSE
                std::clog << "-- Setting empty IRF..." << std::endl;
                std::clog << "-- IRF length: " << _data.size() << std::endl;
#endif
                _irf.reserve(_data.size());
                for (int i = 0; i < _data.size(); i++) _irf.emplace_back(0.0);
                _irf[0] = 1.0;
            } else{
                set_irf(instrument_response_function.data(), instrument_response_function.size());
            }
        }

        // set weights
        if (!weights.empty()) {
            _weights = weights;
        } else {
#if VERBOSE
            std::clog << "-- Assuming Poisson noise" << std::endl;
#endif
            _weights.resize(_data.size());
            for (int i = 0; i < _data.size(); i++)
                _weights[i] = (_data[i] <= 0.0) ? 0.0 : 1. / std::sqrt(_data[i]);
        }
        // set squared weights
        for(int i=0; i<weights.size(); i++){
            _sq_weights.emplace_back(_weights[i] * _weights[i]);
        }

        // set model function
        _model_function.resize(_data.size());

        // sanity check
        if (
                (_data.size() != _time_axis.size()) ||
                (_data.size() != _weights.size()) ||
                (_data.size() != _irf.size())
                ) {
            std::clog << "WARNING: The size of the data, time, weight array, or "
                         "irf do not match" << std::endl;
        }
    }

    /*!
     * Shift an input array by a floating number.
     *
     * @param input[in] the input array
     * @param n_input[in] length of the input array
     * @param output output array
     * @param n_output length of the output array
     * @param shift[in] the shift of the output
     */
    static void shift_array(
            double *input, int n_input,
            double **output_view, int *n_output,
            double shift,
            bool set_outside = true,
            double outside_value = 0.0
    );

    /*!
     * Computes the sum of two arrays considering their respective
     * areal fraction.
     *
     * A weighted sum of two curves is computed. The weighted sum is
     * computed by the area of the first curve and the areal fraction
     * of the second curve. The area of the computed curve equals to
     * the area of the first input curve while the areal fraction of
     * the second input curve will be equal to the specified value in
     * the range specified by the input parameters.
     *
     * @param output the computed output curve (array)
     * @param n_output the number of points in the output
     * @param curve1[in] the first input curve / array
     * @param n_curve1[in] number of points in the first array
     * @param curve2[in] second curve / array
     * @param n_curve2[in] number of points in the second array
     * @param areal_fraction_curve2[in] areal fraction of the second curve in the
     * output array
     * @param start[in] start index used for the area calculation
     * @param stop[in] stop index used for the area calculation
     */
    static void add_curve(
            double **output_view, int *n_output,
            double *curve1, int n_curve1,
            double *curve2, int n_curve2,
            double areal_fraction_curve2,
            int start = 0,
            int stop = -1
    );


    void get_weighted_residuals(double **output_view, int *n_output) {
#if VERBOSE
        std::clog << "Compute weighted residuals..." << std::endl;
        std::clog << "-- points in model function:" << _model_function.size() << std::endl;
        std::clog << "-- points in weights:" << _weights.size() << std::endl;
        std::clog << "-- points in data:" << _data.size() << std::endl;
#endif
        evaluate();
        *n_output = _weighted_residuals.size();
        *output_view = _weighted_residuals.data();
    }

    void evaluate(
            std::vector<double> lifetime_spectrum = std::vector<double>()
    ) {
#if VERBOSE
        std::cout << "evaluate..." << std::endl;
#endif
        if (!lifetime_spectrum.empty()) {
            set_lifetime_spectrum(
                    lifetime_spectrum.data(),
                    lifetime_spectrum.size()
            );
        }
        std::vector<double> lt = lifetime_spectrum;
        if(abs_lifetime_spectrum) {
            for(auto &l: lt){
                l = std::abs(l);
            }
        }
        if (!_is_valid) {
            compute_decay(
                    _model_function.data(),_model_function.size(),
                    _data.data(),_data.size(),
                    _sq_weights.data(),_sq_weights.size(),
                    _time_axis.data(),_time_axis.size(),
                    _irf.data(),_irf.size(),
                    lt.data(),lt.size(),
                    _convolution_start, _convolution_stop,
                    _irf_background_counts, _irf_shift_channels,
                    _areal_scatter_fraction,
                    _period,
                    _constant_background,
                    _total_area,
                    _use_amplitude_threshold,
                    _amplitude_threshold,
                    _correct_pile_up
            );
            if ((_model_function.size() == _data.size()) &&
                (_weights.size() == _data.size())) {
                _weighted_residuals.resize(_data.size());
                for (int i = 0; i < _data.size(); i++) {
                    _weighted_residuals[i] = (_data[i] - _model_function[i]) * _weights[i];
                }
            }
            _is_valid = true;
        }
    }

    /*!
     * Computes the chi2 for the model and the data
     *
     * The "normal" chi2 is the sum of the squared weighted deviations between
     * the data and the model.
     *
     * @param x_min minimum index number of data / model used to compute the chi2
     * @param x_max maximum index number of data / model used to compute the chi2
     * @param type is either neyman or poisson for large count and low count data,
     * respectively.
     * @return the chi2 value
     */
    double get_chi2(
            int x_min = -1,
            int x_max = -1,
            std::string type="poisson"
    ){
#if VERBOSE
        std::cout << "CHI2" << std::endl;
#endif
        if(x_min < 0){
            x_min = std::max(this->x_min, 0);
        } else{
            x_min = std::max(x_min, 0);
        }
        if(x_max < 0){
            if(this->x_max < 0){
                x_max = (int)_data.size();
            } else{
                x_max = std::min(this->x_max, (int)_data.size());
            }
        }
        if(!is_valid()) {
            evaluate();
#if VERBOSE
            std::cout << "-- evaluate" << std::endl;
#endif
        }
        double v = statistics::chi2_counting(
                _data,
                _model_function,
                x_min, x_max, type
        );
#if VERBOSE
        std::cout << "-- x_min: " << x_min << std::endl;
        std::cout << "-- x_max: " << x_max << std::endl;
        std::cout << "-- type: " << type << std::endl;
        std::cout << "-- chi2: " << v << std::endl;
#endif
        return v;
    }

    /*!
     * Convenience method to update parameters that are frequently changed.
     *
     * @param irf_background
     * @param irf_shift_channels
     * @param areal_scatter_fraction
     * @param constant_background
     * @param lifetime_spectrum
     * @param n_lifetime_spectrum
     */
    void set(
        double irf_background = 0.0,
        double irf_shift_channels = 0.0,
        double areal_scatter_fraction = 0.0,
        double constant_background = 0.0,
        double *lifetime_spectrum = nullptr, int n_lifetime_spectrum = -1
    ){
        set_irf_background_counts(irf_background);
        set_irf_shift_channels(irf_shift_channels);
        set_constant_background(constant_background);
        set_lifetime_spectrum(lifetime_spectrum, n_lifetime_spectrum);
        set_areal_scatter_fraction(areal_scatter_fraction);
    }

};


#endif //TTTRLIB_DECAY_H
