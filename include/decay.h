#ifndef FIT2X_DECAY_H
#define FIT2X_DECAY_H

#include <cmath>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <numeric>
#include <algorithm>
#include <memory> /* shared_ptr */
#include "omp.h"

#include "fsconv.h"
#include "tttrlib/tttr.h"
#include "statistics.h"



class Decay {

private:

    /// Used to keep track of irf shifts and background
    /// if false get_corrected_irf will recompute the irf
    bool _irf_is_corrected = true;

    /// The model function is multiplied by this vector is _use_linearization is true
    std::vector<double> _linearization_table;

    /// If set to true multiply the linearization to the model function
    bool _use_linearization = false;

    /// If set to true scales the model to the data (default false)
    bool _scale_model_to_data;

    /// If set to true uses the background corrected IRF as scatter
    bool _use_corrected_irf_as_scatter;

    /// Acquisition time of the experiment in seconds
    double _acquisition_time;

    /// Dead time of the instrument in units of the lifetime (usually nanoseconds)
    double _instrument_dead_time;

    /// Stores the weighted residuals
    std::vector<double> _weighted_residuals = std::vector<double>();

    /// Proportional to the area of the model. If negative, the area is scaled
    /// to the data.
    double _number_of_photons = -1;

    /// The constant offset in the model
    double _constant_offset = 0.0;

    /// The fraction of the irf (scattered light) in the model decay
    double _scatter_fraction = 0.0;

    /// The time shift of the irf
    double _irf_shift_channels = 0.0;

    /// Background counts of irf. This number is sutracted from the counts
    /// in each channel before convolution.
    double _irf_background_counts = 0.0;

    /// The repetiotion period (usually in nano seconds)
    double _excitation_period = 100.;

    /// The background corrected instrument response function
    std::vector<double> _corrected_irf = std::vector<double>();

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
    bool _use_pile_up_correction = false;

    /// Convolution range
    int _convolution_start = 0;
    int _convolution_stop = -1;

    /// Amplitude threshold used to discriminate low populated lifetimes
    double _amplitude_threshold = 1e10;

    /// If set to true discriminates low populated lifetimes in convolution.
    bool _use_amplitude_threshold = false;

    /// If set to true take abs values of lifetime spectrum
    bool _abs_lifetime_spectrum = true;

    /// The range on which the scoring function is evaluated
    int _score_range_min = -1;
    int _score_range_max = -1;

    /// The method used for convolution
    int _convolution_method = 0;

    /// Used to set if the decay is 'valid'. A decay is invalid if the
    /// computed model, wres, etc. does not correspond to the input parameters.
    /// It should not be decided by outside if the decay is valid.
    void set_is_valid(bool v){
        _is_valid = v;
    }

protected:

    static void scale_model(
            bool scale_model_to_data,
            double number_of_photons,
            int convolution_start, int convolution_stop,
            double constant_background,
            double* model,
            double* data,
            double* squared_data_weights
    );

public:

    std::vector<int> get_score_range(){
        return std::vector<int>({_score_range_min, _score_range_max});
    }

    void set_score_range(int min, int max){
        _score_range_min = min;
        _score_range_max = max;
    }


    bool get_is_valid() const {
        return _is_valid;
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
        set_is_valid(false);
        _convolution_method = v;
    }

    int get_convolution_method(){
        return _convolution_method;
    }

    void set_abs_lifetime_spectrum(bool v){
        set_is_valid(false);
        _abs_lifetime_spectrum = v;
    }

    bool get_abs_lifetime_spectrum() const{
        return _abs_lifetime_spectrum;
    }

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
     * @param squared_data_weights[in] The squared weights of the data points. The data
     * weights are used to scale the model function to the data (usually the data
     * weights is Poissonian)
     * @param n_data_weights[in] The number of weights
     * @param time_axis[in] The time axiss for which the fluorescence decay is computed
     * for
     * @param n_time_axis[in] The number of points in the time axis
     * @param irf_histogram[in] The instrument response function
     * used for convolution
     * @param n_irf_histogram[in] The number of points in the
     * instrument response function
     * @param lifetime_spectrum[in] The lifetime spectrum. A lifetime specturm is an
     * interleaved array of amplitudes and fluorescence lifetimes, .e.g,
     * (amplitude1, lifetime1, amplitude2, lifetime2, ...)
     * @param n_lifetime_spectrum[in] The number of points in the fluorescence lifetime
     * spectrum
     * @param convolution_start[in] The start of the convolution
     * @param convolution_stop[in] The stop of the convoution. The decay will be computed in the
     * range [start, stop]
     * @param irf_background_counts[in] The background counts of the instrument
     * response function. This number will be subtracted from the IRF before
     * convolution.
     * @param irf_shift_channels[in] The number of micro time channels the IRF will be
     * shifted before the fluorescence lifetimes are convoluted with the IRF.
     * @param scatter_fraction[in] The fraction (integrated fraction), i.e.,
     * the area the scattered light will have in the computed decay.
     * @param excitation_period[in] The excitation period (in units of the fluorescence
     * lifetime, usually nanoseconds) that was used to excite the sample
     * @param constant_offset A constant offset that is added to the fluorescence
     * decay.
     * @param number_of_photons[in] the area to which the model fluorescence decay is scaled
     * to. If this value is negative (default value is -1) the model fluorescence
     * decay is scaled to the data in the range defined by the parameters start,
     * stop
     * @param use_amplitude_threshold[in] if set to true (default is false) a
     * threshold will be applied. If this parameter is set to true, fluorescence
     * lifetimes with a amplitude that is smaller than a threshold will be not
     * considered.
     * @param amplitude_threshold The threshold that is used to discriminate
     * fluorescence lifetimes that are smaller.
     * @param use_pile_up_correction if set to true (default is false) pile up will be added
     * to the model function.
     * @param instrument_dead_time the dead time of the instrument (used for
     * pile up (in units of the lifetime, usually nano seconds)
     * @param acquisition_time the total time the acquisition of the decay in
     * seconds.
     * @param add_corrected_irf_as_scatter if set to true (default is false) the background
     * corrected irf will be added as scatter to the decay. If this is false the
     * irf prior to a background and shift corrected irf is added as a scatter
     * fraction.
     * @param scale_model_to_data if set to true (default is true) the model is
     * either scaled to the area s provided by total_area (if total_area is larger
     * then zero) or to the total number of experimental counts.
     */
    static int compute_decay(
            double* model_function, int n_model_function,
            double* data, int n_data,
            double* squared_data_weights, int n_squared_data_weights,
            double* time_axis, int n_time_axis,
            double* irf_histogram, int n_irf_histogram,
            double* lifetime_spectrum, int n_lifetime_spectrum,
            double* scatter = nullptr, int n_scatter = -1,
            int convolution_start = 0, int convolution_stop = -1,
            double scatter_fraction = 0.0,
            double excitation_period = 1000.0,
            double constant_offset = 0.0,
            double number_of_photons=1,
            bool use_amplitude_threshold = false,
            double amplitude_threshold = 1e10,
            bool use_pile_up_correction = false,
            double instrument_dead_time=120.0,
            double acquisition_time=1e9,
            bool add_corrected_irf_as_scatter = false,
            bool scale_model_to_data = false,
            int convolution_method = 0,
            bool take_abs_of_lifetime_spectrum = false,
            bool use_linearization = false,
            double* linearization= nullptr, int n_linearization = -1
    );


    void set_weights_by_data(std::vector<double> data){
#if VERBOSE_FIT2X
        std::clog << "-- Using Poisson noise" << std::endl;
#endif
        set_is_valid(false);
        _sq_weights.resize(data.size());
        _weights.resize(data.size());
        for (size_t i = 0; i < data.size(); i++){
            double w = (data[i] <= 0.0) ? 0.0 : 1. / std::sqrt(data[i]);
            _weights[i] = w;
            _sq_weights[i] = w * w;
        }
    }

    void set_data(double *input, int n_input) {
        set_is_valid(false);
        _model_function.resize(n_input);
        _data.resize(n_input);
        _data.assign(input, input + n_input);
        set_weights_by_data(_data);
    }

    void get_data(double **output_view, int *n_output) {
        *output_view = _data.data();
        *n_output = _data.size();
    }

    void set_use_amplitude_threshold(bool v) {
        set_is_valid(false);
        _use_amplitude_threshold = v;
    }

    bool get_use_amplitude_threshold() const {
        return _use_amplitude_threshold;
    }

    void set_amplitude_threshold(double v) {
        set_is_valid(false);
        _amplitude_threshold = v;
    }

    double get_amplitude_threshold() const {
        return _amplitude_threshold;
    }

    void set_number_of_photons(double v) {
        set_is_valid(false);
        _number_of_photons = v;
    }

    double get_number_of_photons(){
        if(_scale_model_to_data){
            double re = 0.0;
            double* data; int n_data;
            get_data(&data, &n_data);
            for(int i=get_convolution_start(); i<get_convolution_stop();i++){
                re += data[i];
            }
            return re;
        }
        return _number_of_photons;
    }

    void set_excitation_period(double v) {
        set_is_valid(false);
        _excitation_period = v;
    }

    double get_excitation_period() const {
        return _excitation_period;
    }

    void set_irf_shift_channels(double v) {
        set_is_valid(false);
        double n_channels = std::max(1., (double)_irf.size());
        _irf_shift_channels = std::fmod(v, n_channels);
        _irf_is_corrected = false;
    }

    double get_irf_shift_channels() const {
        return _irf_shift_channels;
    }

    void set_scatter_fraction(double v) {
        set_is_valid(false);
        _scatter_fraction = v; //std::min(1., std::max(0.0, v));
    }

    double get_scatter_fraction() const {
        return _scatter_fraction;
    }

    void set_constant_offset(double v) {
        set_is_valid(false);
        _constant_offset = v;
    }

    double get_constant_offset() const {
        return _constant_offset;
    }

    void set_convolution_start(int v) {
        set_is_valid(false);
        _convolution_start = std::max(0, v);
    }

    int get_convolution_start() const {
        return std::max(0, _convolution_start);
    }

    void set_convolution_stop(int v) {
        set_is_valid(false);
        _convolution_stop = v;
    }

    int get_convolution_stop() const {
        int v = (_convolution_stop < 0) ? (int) _data.size() : std::min((int) _data.size(), _convolution_stop);
        return v;
    }

    void set_convolution_range(std::vector<int> v) {
        // set_is_valid(false); // this is already set by set_convolution_start, set_convolution_stop
        set_convolution_start(v[0]);
        set_convolution_stop(v[1]);
    }

    std::vector<int> get_convolution_range() const {
        return std::vector<int>({_convolution_start, _convolution_stop});
    }

    void set_use_pile_up_correction(bool v) {
        set_is_valid(false);
        _use_pile_up_correction = v;
    }

    bool get_use_pile_up_correction() const {
        return _use_pile_up_correction;
    }

    void set_irf(double *input, int n_input) {
#if VERBOSE_FIT2X
        std::clog << "-- set_irf" << std::endl;
#endif
        set_is_valid(false);
        if(n_input>0){
            _irf.resize(n_input);
            _irf.assign(input, input + n_input);
            _irf_is_corrected = false;
        }
        else {
            if (!_data.empty()) {
#if VERBOSE_FIT2X
                std::clog << "-- Setting empty IRF..." << std::endl;
                std::clog << "-- IRF length: " << _data.size() << std::endl;
#endif
                _irf.resize(_data.size());
                std::fill(_irf.begin(), _irf.end(), 0.0);
                _irf[0] = 1.0;
            }
        }
    }

    void get_irf(double **output_view, int *n_output) {
        *output_view = _irf.data();
        *n_output = _irf.size();
    }

    void set_linearization(double *input, int n_input) {
#if VERBOSE_FIT2X
        std::clog << "-- set_linearization" << std::endl;
#endif
        set_is_valid(false);
        _linearization_table.resize(n_input);
        _linearization_table.assign(input, input + n_input);
        // Sanity check
        if(n_input < (int) _data.size()){
            std::clog << "-- WARNING: linearization too short filling with ones." << std::endl;
            while(_linearization_table.size() < _data.size()){
                _linearization_table.emplace_back(1.0);
            }
        }
    }

    void get_linearization(double **output_view, int *n_output) {
        *output_view = _linearization_table.data();
        *n_output = _linearization_table.size();
    }

    bool get_use_linearization() const{
        return _use_linearization;
    }

    void set_use_linearization(bool v){
        set_is_valid(false);
        _use_linearization = v;
    }

    void get_corrected_irf(double **output_view, int *n_output) {
        if (!_irf_is_corrected) {
            // correct irf for background counts
            _corrected_irf.resize(_irf.size());
            for (size_t i = 0; i < _irf.size(); i++)
                _corrected_irf[i] = std::max(_irf[i] - _irf_background_counts, 0.0);
            // shift irf
            _corrected_irf = shift_array(
                    _corrected_irf.data(),
                    _corrected_irf.size(),
                    _irf_shift_channels
            );
            _irf_is_corrected = true;
        }
        *output_view = _corrected_irf.data();
        *n_output = _corrected_irf.size();
    }

    void get_model(double **output_view, int *n_output){
        if (!_is_valid) evaluate();
        *output_view = _model_function.data();
        *n_output = _model_function.size();
    }

    void set_lifetime_spectrum(double *input, int n_input = -1) {
        set_is_valid(false);
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

    void set_data_weights(double *input, int n_input) {
        set_is_valid(false);
#if VERBOSE_FIT2X
        std::clog << "-- Setting weights..." << std::endl;
#endif
        if (n_input > 0) {
            _weights.resize(n_input);
            _weights.assign(input, input+n_input);
            _sq_weights.clear();
            for(auto &v: _weights) _sq_weights.emplace_back(v * v);
        } else {
#if VERBOSE_FIT2X
            std::clog << "-- WARNING: No weights provided" << std::endl;
#endif
            set_weights_by_data(_data);
        }
    }

    void get_data_weights(double **output_view, int *n_output) {
        *output_view = _weights.data();
        *n_output = _weights.size();
    }

    void set_time_axis(double *input, int n_input) {
        set_is_valid(false);
        _time_axis.resize(n_input);
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
        _irf_is_corrected = false;
        set_is_valid(false);
    }

    double get_irf_background_counts() const {
        return _irf_background_counts;
    }

    void set_instrument_dead_time(
            double instrument_dead_time
    ) {
        set_is_valid(false);
        _instrument_dead_time = instrument_dead_time;
    }

    double get_instrument_dead_time() const {
        return _instrument_dead_time;
    }

    void set_acquisition_time(
            double acquisition_time
    ) {
        set_is_valid(false);
        if(acquisition_time < 0) acquisition_time = 1e12;
        _acquisition_time = acquisition_time;
    }

    double get_acquisition_time() const {
        return _acquisition_time;
    }

    void set_use_corrected_irf_as_scatter(bool v) {
        set_is_valid(false);
        _use_corrected_irf_as_scatter = v;
    }

    bool get_use_corrected_irf_as_scatter() const {
        return _use_corrected_irf_as_scatter;
    }

    void set_scale_model_to_data(bool v) {
        set_is_valid(false);
        _scale_model_to_data = v;
    }

    bool get_scale_model_to_data() const {
        return _scale_model_to_data;
    }

    void set_tttr_data(std::shared_ptr<TTTR> tttr_data, int tttr_micro_time_coarsening){
        set_is_valid(false);
        double *hist; int n_hist;
        double *time; int n_time;
        TTTR::compute_microtime_histogram(
                tttr_data.get(),
                &hist, &n_hist,
                &time, &n_time,
                tttr_micro_time_coarsening
        );
        set_data(hist, n_hist);
        set_time_axis(time, n_time);
        free(hist); free(time);
        _time_axis.clear();
        _time_axis.reserve(_data.size());
        double micro_time_resolution = tttr_data->get_header().micro_time_resolution / tttr_micro_time_coarsening;
        for(size_t i = 0; i < _data.size(); i++) _time_axis.emplace_back(i * micro_time_resolution);
    }

    void set_tttr_irf(std::shared_ptr<TTTR> tttr_irf, int tttr_micro_time_coarsening){
#if VERBOSE_FIT2X
        std::clog << "-- Setting IRF from TTTR..." << std::endl;
#endif
        set_is_valid(false);
        double *hist; int n_hist;
        double *time; int n_time;
        TTTR::compute_microtime_histogram(
                tttr_irf.get(),
                &hist, &n_hist,
                &time, &n_time,
                tttr_micro_time_coarsening
        );
        set_irf(hist, n_hist);
        free(time); free(hist);
    }

    /*!
     *
     * @param tttr_data pointer to TTTR object that is used to construct a decay
     * histogram
     * @param tttr_micro_time_coarsening an (optional) integer by which the micro times
     * are divided to coarsen the time axis (default is 1)
     * @param data the data to which the decay is fitted
     * @param time_axis the time axis that belongs to the data
     * @param data_weights the weights of the data points. If the weights are not provided
     * (nullptr / None) the weights are computed assuming Poisson noise.
     * @param irf_histogram The instrument response function (IRF)
     * that is used for convolution. If no IRF is provided
     * @param convolution_start The start index in the IRF used for convolution. Points in the
     * IRF before the start index are not used for convolution.
     * @param convolution_stop The stop index in the IRF used for convolution. Points beyond the
     * stop index are not convolved.
     * @param use_amplitude_threshold If this is set to true (default value is true)
     * the values that are smaller then a specified threshold are omitted
     * @param amplitude_threshold The amplitude threshold that is used if the
     * parameter use_amplitude_threshold is set to true (the default value is 1e10)
     * @param use_pile_up_correction If this is set to true (the default value is false)
     * the convolved model function is 'piled up' to match pile up artifacts in the
     * data.
     * @param excitation_period the repetition period, .i.e, the time between subsequent
     * excitation pulses.
     */
    Decay(
            std::vector<double> data = std::vector<double>(),
            std::vector<double> time_axis = std::vector<double>(),
            std::vector<double> irf_histogram = std::vector<double>(),
            double constant_offset = 0.0,
            double irf_background_counts = 0.0,
            double scatter_fraction = 0.0,
            double irf_shift_channels=0.0,
            std::vector<double> lifetime_spectrum = std::vector<double>(),
            std::vector<int>  convolution_range = std::vector<int>({0, -1}),
            std::vector<int> score_range = std::vector<int>({0, -1}),
            double excitation_period = 100.0,
            bool scale_model_to_data = true,
            double number_of_photons=1,
            std::vector<double> linearization = std::vector<double>(),
            bool use_linearization = false,
            bool use_pile_up_correction = false,
            bool use_corrected_irf_as_scatter = false,
            bool use_amplitude_threshold = false,
            double amplitude_threshold = 1e-9,
            double acquisition_time = 1e6,
            double instrument_dead_time = 1e-9,
            int convolution_method = 0,
            std::vector<double> data_weights = std::vector<double>(),
            std::shared_ptr<TTTR> tttr_data = nullptr,
            std::shared_ptr<TTTR> tttr_irf = nullptr,
            int tttr_micro_time_coarsening = 1
    ) {
#if VERBOSE_FIT2X
        std::clog << "NEW DECAY" << std::endl;
#endif
        set(
                data,
                time_axis,
                irf_histogram,
                constant_offset,
                irf_background_counts,
                scatter_fraction,
                irf_shift_channels,
                lifetime_spectrum,
                convolution_range,
                score_range,
                excitation_period,
                scale_model_to_data,
                number_of_photons,
                linearization,
                use_linearization,
                use_pile_up_correction,
                use_corrected_irf_as_scatter,
                use_amplitude_threshold,
                amplitude_threshold,
                acquisition_time,
                instrument_dead_time,
                convolution_method,
                data_weights,
                tttr_data,
                tttr_irf,
                tttr_micro_time_coarsening,
                true, // to make sure that the length of the irf matches the data
                true // to make sure that the time axis is filled
        );
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
    static std::vector<double> shift_array(
            double *input, int n_input,
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
     * modifies curve1 inplace
     *
     * @param output the computed output curve (array)
     * @param n_output the number of points in the output
     * @param curve1[in, out] the first input curve / array
     * @param n_curve1[in, out] number of points in the first array
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

    static double compute_score(
            double *data, int n_data,
            double *time_axis = nullptr, int n_time_axis = -1,
            double *irf_histogram = nullptr, int n_irf_histogram = -1,
            double constant_offset = 0.0,
            double irf_background_counts = 0.0, double scatter_fraction = 0.0,
            double irf_shift_channels = 0.0, double *lifetime_spectrum = nullptr,
            int n_lifetime_spectrum = -1,
            std::vector<int> convolution_range = std::vector<int>({0, -1}),
            std::vector<int> score_range = std::vector<int>({0, -1}),
            double excitation_period = 100.0,
            bool scale_model_to_data = true, double number_of_photons = 1,
            double *linearization = nullptr, int n_linearization = -1,
            bool use_linearization = false,
            bool use_pile_up_correction = false, bool use_corrected_irf_as_scatter = false,
            bool use_amplitude_threshold = false, double amplitude_threshold = 1e10,
            double acquisition_time = 1e9, double instrument_dead_time = 1e-9,
            int convolution_method = 0,
            double *data_weights = nullptr, int n_data_weights = -1,
            std::string score_type = "poisson",
            bool take_abs_of_lifetime_spectrum = true
    ) {
        int n_model = std::min({n_irf_histogram, n_time_axis, n_data});
        std::vector<double> model_function(n_model, 0.0);
        std::vector<double> irf(n_irf_histogram, 0.0);

        std::vector<double> sq_weights(n_data_weights, 1.0); // by default one to avoid div by zero
        for(int i=0; i < n_data_weights; i++) sq_weights[i] = data_weights[i] * data_weights[i];

        for(int i=0; i < n_irf_histogram; i++) irf[i] = std::max(0.0, irf_histogram[i] - irf_background_counts);
        irf = shift_array(irf.data(), irf.size(), irf_shift_channels);
        // returns 0 if success
        int r = compute_decay(
                model_function.data(),
                model_function.size(),
                data, n_data,
                sq_weights.data(), sq_weights.size(),
                time_axis, n_time_axis,
                irf.data(), irf.size(),
                lifetime_spectrum, n_lifetime_spectrum,
                irf_histogram, n_irf_histogram,
                convolution_range[0], convolution_range[1],
                scatter_fraction,
                excitation_period,
                constant_offset,
                number_of_photons,
                use_amplitude_threshold, amplitude_threshold,
                use_pile_up_correction,
                instrument_dead_time, acquisition_time,
                use_corrected_irf_as_scatter,
                scale_model_to_data,
                convolution_method,
                take_abs_of_lifetime_spectrum,
                use_linearization,
                linearization, n_linearization
        );
        if(r >= 0){
            auto decay_histogram = std::vector<double>(data, data+n_data);
            double score = statistics::chi2_counting(
                    decay_histogram,
                    model_function,
                    score_range[0], score_range[1],
                    score_type
            );
            return score;
        } else{
            return INFINITY;
        }
    }

    void get_weighted_residuals(double **output_view, int *n_output) {
        double *m; int nm; get_model(&m, &nm);
        double *d; int nd; get_data(&d, &nd);
        double *w; int nw; get_data_weights(&w, &nw);
#if VERBOSE_FIT2X
        std::clog << "Compute weighted residuals..." << std::endl;
        std::clog << "-- points in model function: " << nm << std::endl;
        std::clog << "-- points in weights: " << nw << std::endl;
        std::clog << "-- points in data: " << nd << std::endl;
#endif
        if (nm == nd && nd == nw) _weighted_residuals.resize(nd);
        for (int i = 0; i < nd; i++)
            _weighted_residuals[i] = (d[i] - m[i]) * w[i];
        *n_output = _weighted_residuals.size();
        *output_view = _weighted_residuals.data();
    }

    void evaluate() {
#if VERBOSE_FIT2X
        std::clog << "evaluate..." << std::endl;
#endif
        if (!get_is_valid()) {
            double* irf; int n_irf;
            get_corrected_irf(&irf, &n_irf);
            double* scatter; int n_scatter;
            get_irf(&scatter, &n_scatter);
            int s = compute_decay(
                    _model_function.data(), _model_function.size(),
                    _data.data(), _data.size(),
                    _sq_weights.data(),
                    _sq_weights.size(),
                    _time_axis.data(),
                    _time_axis.size(),
                    irf, n_irf,
                    _lifetime_spectrum.data(), _lifetime_spectrum.size(),
                    scatter, n_scatter,
                    get_convolution_start(), get_convolution_stop(),
                    get_scatter_fraction(),
                    get_excitation_period(),
                    _constant_offset,
                    _number_of_photons, // the getter calls this methods
                    get_use_amplitude_threshold(),
                    get_amplitude_threshold(),
                    get_use_pile_up_correction(),
                    get_instrument_dead_time(),
                    _acquisition_time,
                    get_use_corrected_irf_as_scatter(),
                    get_scale_model_to_data(),
                    get_convolution_method(),
                    get_abs_lifetime_spectrum(),
                    get_use_linearization(),
                    _linearization_table.data(), _linearization_table.size()
            );
            set_is_valid((s >= 0));
        }
    }

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
        double m0_irf = std::accumulate(
                irf_histogram.begin(),
                irf_histogram.end(),
                0.0);
        double m1_irf = 0.0;
        for(size_t i = 0; i<irf_histogram.size(); i++)
            m1_irf += i * irf_histogram[i];

        double mu0 = std::accumulate(
                decay_histogram.begin(),
                decay_histogram.end(),
                0.0);
        double mu1 = 0.0;
        for(size_t i = 0; i<decay_histogram.size(); i++)
            mu1 += i * decay_histogram[i];

        double g1 = mu0 / m0_irf;
        double g2 = (mu1 - g1 * m1_irf) / m0_irf;
        double tau1 = g2 / g1 * micro_time_resolution;
        return tau1;
    }

    /*!
     *
     * @return mean lifetime compute by the methods of momements (first moment)
     */
    double get_mean_lifetime(){
        double dt = 1.0;
        if(_time_axis.size() > 1) dt = _time_axis[1] - _time_axis[0];
        return compute_mean_lifetime(_irf, _data, dt);
    }

    /*!
     * Computes the chi2 for the model and the data
     *
     * The "normal" chi2 is the sum of the squared weighted deviations between
     * the data and the model.
     *
     * @param x_min minimum index number of data / model used to compute the chi2
     * @param x_max maximum index number of data / model used to compute the chi2
     * @param score_type is either neyman or poisson for large count and low count data,
     * respectively, pearson, gauss, cnp
     * @return the chi2 value
     */
    double get_score(
            int x_min = -1,
            int x_max = -1,
            std::string score_type= "poisson"
    ){
#if VERBOSE_FIT2X
        std::clog << "CHI2" << std::endl;
        std::clog << "-- data range: " << x_min << ", " << x_max << std::endl;
#endif
        x_min = (x_min < 0) ? std::max(this->_score_range_min, 0) : x_min;
        x_max = (x_max < 0) ? (int)_data.size(): std::min(x_max, (int)_data.size());
        double v = 0.0;
        if(score_type != "normal") {
            double* m; int nm; get_model(&m, &nm);
            double* d; int nd; get_data(&d, &nd);
            auto data = std::vector<double>(d, d+nd);
            auto model = std::vector<double>(m, m+nm);
            v = statistics::chi2_counting(data, model, x_min, x_max, score_type);
        } else{
            double* wres; int nwres;
            get_weighted_residuals(&wres, &nwres);
            auto sr = get_score_range();
            for(int i=sr[0];i<sr[1];i++) v += wres[i] * wres[i];
        }
#if VERBOSE_FIT2X
        std::clog << "-- x_min: " << x_min << std::endl;
        std::clog << "-- x_max: " << x_max << std::endl;
        std::clog << "-- chi2: " << v << std::endl;
#endif
        if(get_is_valid()) return v;
        else return INFINITY;
    }

    /*!
     * Update parameters
     */
    void set(
            std::vector<double> data = std::vector<double>(),
            std::vector<double> time_axis = std::vector<double>(),
            std::vector<double> irf_histogram = std::vector<double>(),
            double constant_offset = 0.0,
            double irf_background_counts = 0.0,
            double scatter_fraction = 0.0,
            double irf_shift_channels = 0.0,
            std::vector<double> lifetime_spectrum = std::vector<double>(),
            std::vector<int> convolution_range = std::vector<int>({0, -1}),
            std::vector<int> score_range = std::vector<int>({0, -1}),
            double excitation_period = 100.0,
            bool scale_model_to_data = true,
            double number_of_photons = 1,
            std::vector<double> linearization = std::vector<double>(),
            bool use_linearization = false,
            bool use_pile_up_correction = false,
            bool use_corrected_irf_as_scatter = false,
            bool use_amplitude_threshold = false,
            double amplitude_threshold = 1e-9,
            double acquisition_time = 1e9,
            double instrument_dead_time = 1e-9,
            int convolution_method = 0,
            std::vector<double> data_weights = std::vector<double>(),
            std::shared_ptr<TTTR> tttr_data = nullptr,
            std::shared_ptr<TTTR> tttr_irf = nullptr,
            int tttr_micro_time_coarsening = 1,
            bool force_fill_irf = false,
            bool force_time_axis = false
    ) {
        // set data
        if (tttr_data != nullptr) {
            excitation_period = (tttr_data == nullptr) ?
                                excitation_period :
                                tttr_data->get_header().macro_time_resolution;
            set_tttr_data(tttr_data, tttr_micro_time_coarsening);
        } else{
            if(!data.empty()) set_data(data.data(), data.size());
            if(!time_axis.empty()) set_time_axis(time_axis.data(), time_axis.size());
            else if(force_time_axis){
                auto t = std::vector<double>(data.size());
                for(uint32_t i=0; i < t.size(); i++) t[i] = i;
                set_time_axis(t.data(), t.size());
            }
        }

        // set irf
        if (tttr_irf != nullptr)
            set_tttr_irf(tttr_irf, tttr_micro_time_coarsening);
        else if(!irf_histogram.empty())
            set_irf(irf_histogram.data(), irf_histogram.size());
        else if(force_fill_irf){
            // This sets the irf to the length of the data. The first element
            // in the irf is a one.
            set_irf(nullptr, 0);
        }

        // set linearization
        if(!linearization.empty()){
            set_linearization(linearization.data(), linearization.size());
            set_use_linearization(use_linearization);
        }
        else{
            set_use_linearization(false);
        }

        // set weights
        if(!data_weights.empty()) set_data_weights(data_weights.data(), data_weights.size());

        if (score_range.empty() < 1) {
            score_range.emplace_back(0);
            score_range.emplace_back(-1);
        } else if (score_range.size() < 2) {
            score_range.emplace_back(-1);
        }
        // convolution method
        set_convolution_method(convolution_method);

#if VERBOSE_FIT2X
        std::clog << "SET" << std::endl;
        std::clog << "irf_background_counts: " << irf_background_counts << std::endl;
        std::clog << "irf_shift_channels: " << irf_shift_channels << std::endl;
        std::clog << "scatter_fraction: " << scatter_fraction << std::endl;
        std::clog << "constant_offset: " << constant_offset << std::endl;
        std::clog << "number_of_photons: " << number_of_photons << std::endl;
        std::clog << "acquisition_time: " << acquisition_time << std::endl;
        std::clog << "instrument_dead_time: " << instrument_dead_time << std::endl;
        std::clog << "convolution_range: " << _convolution_start << ", " << _convolution_stop << std::endl;
        std::clog << "excitation_period: " << excitation_period << std::endl;
        std::clog << "scale_model_to_data: " << scale_model_to_data << std::endl;
        std::clog << "score_range: " << score_range[0] << ", " << score_range[1] << std::endl;
        std::clog << "use_corrected_irf_as_scatter: " << use_corrected_irf_as_scatter << std::endl;
        std::clog << "amplitude_threshold: " << amplitude_threshold << std::endl;
        std::clog << "use_amplitude_threshold: " << use_amplitude_threshold << std::endl;
        std::clog << "use_linearization: " << use_linearization << std::endl;
        std::clog << "---" << std::endl;
#endif
        set_number_of_photons(number_of_photons);
        set_irf_background_counts(irf_background_counts);
        set_irf_shift_channels(irf_shift_channels);
        set_constant_offset(constant_offset);
        if(!lifetime_spectrum.empty())
            set_lifetime_spectrum(lifetime_spectrum.data(), lifetime_spectrum.size());
        set_scatter_fraction(scatter_fraction);
        set_acquisition_time(acquisition_time);
        set_instrument_dead_time(instrument_dead_time);
        set_convolution_range(convolution_range);
        set_excitation_period(excitation_period);
        set_scale_model_to_data(scale_model_to_data);
        set_score_range(score_range[0], score_range[1]);
        set_use_corrected_irf_as_scatter(use_corrected_irf_as_scatter);
        set_amplitude_threshold(amplitude_threshold);
        set_use_amplitude_threshold(use_amplitude_threshold);
        set_use_pile_up_correction(use_pile_up_correction);
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

};


#endif //FIT2X_DECAY_H

