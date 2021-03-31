#ifndef FIT2X_DECAY_H
#define FIT2X_DECAY_H

#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <cstring> /* strcmp */
#include <memory> /* shared_ptr */
#include <limits> /* std::numeric_limits */

#include "omp.h"
#include "tttrlib/TTTR.h"

#include "fsconv.h"
#include "statistics.h"



class Decay {

private:

    /// The size of the decay (number of channels)
    int _size = 0;

    /// Used to keep track if corrected irf needs an update
    bool _irf_is_corrected = true;

    /// The model function is multiplied by this vector is _use_linearization is true
    std::vector<double> _linearization_table;

    /// If set to true multiply the linearization to the model function
    bool _use_linearization = false;

    /// If set to true scales the model to the data (default false)
    bool _scale_model_to_data = true;

    /// If set to true uses the background corrected IRF as scatter
    bool _use_corrected_irf_as_scatter = true;

    /// Acquisition time of the experiment in seconds
    double _acquisition_time = std::numeric_limits<double>::max();

    /// Dead time of the instrument in units of the lifetime (usually nanoseconds)
    double _instrument_dead_time = std::numeric_limits<double>::epsilon();

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

    /// The excitation repetition period (usually in nano seconds)
    double _excitation_period = std::numeric_limits<double>::max();

    /// The instrument response function
    std::vector<double> _irf = std::vector<double>();

    /// The background and shift corrected irf
    std::vector<double> _corrected_irf = std::vector<double>();

    /// The experimental histogram
    std::vector<double> _data = std::vector<double>();

    /// The weights of the experimental data
    std::vector<double> _weights = std::vector<double>();

    /// The squared weights of the experimental data
    std::vector<double> _sq_weights = std::vector<double>();

    /// The time axis of the data and model
    std::vector<double> _time_axis = std::vector<double>();

    /// model function
    std::vector<double> _model_function = std::vector<double>();

    /// The lifetime spectrum
    std::vector<double> _lifetime_spectrum = std::vector<double>();

    /// false if model function needs to be updated
    bool _is_valid = false;

    /// If true pile up distortion will added to model function
    bool _use_pile_up_correction = false;

    /// Convolution range
    int _convolution_start = 0;
    int _convolution_stop = 0;

    /// Threshold used to discriminate lifetimes with small amplitudes
    double _amplitude_threshold = std::numeric_limits<double>::epsilon();

    /// If true lifetimes with small amplitudes are discriminated
    bool _use_amplitude_threshold = false;

    /// If true absolute values of lifetime spectrum is used to compute model function
    bool _abs_lifetime_spectrum = true;

    /// The range on which the scoring function is evaluated
    int _score_range_min = 0;
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

    /*!
     * Scales the model in the specified range
     *
     * The model is either scaled to the data or to a specified number of photons
     *
     * @param scale_model_to_data[in] If set to true the model is scaled to the
     * data
     * @param number_of_photons[in] If scale_model_to_data is fales the model is
     * scaled to the specified number of photons
     * @param start[in] Specifies the start index of the model range
     * @param stop[in] Specifies the stop index of the model range
     * @param constant_background[in] number of constant background photons in
     * each histogram bin.
     * @param model[in, out] Model function that is scaled. The model function is
     * changed inplace.
     * @param data[in] Data array to which the model function can be scaled
     * @param squared_data_weights[in] squared weights of the data
     */
    static void scale_model(
            bool scale_model_to_data,
            double number_of_photons,
            int start, int stop,
            double constant_background,
            double* model,
            double* data,
            double* squared_data_weights
    );

    /// Equalize the length of the
    void equalize_length(){
        size_t n_min = std::max({
            _irf.size(),
            _model_function.size(),
            _data.size(),
            _linearization_table.size()
        });
        resize(n_min);
#if VERBOSE_FIT2X
        std::clog << "-- equalize_length: " << n_min << std::endl;
#endif
    }

    void set_weights_by_data(const char* noise_model = "poisson"){
#if VERBOSE_FIT2X
        std::clog << "-- Using Poisson noise" << std::endl;
#endif
        set_is_valid(false);
        if(strcmp(noise_model, "poisson") == 0){
            for (int i = 0; i < size(); i++){
                double w = (_data[i] <= 0.0) ? 0.0 : 1. / std::sqrt(_data[i]);
                _weights[i] = w;
                _sq_weights[i] = w * w;
            }
        }
    }

public:

    double get_time_resolution(){
        double dt = 1.0;
        if(_time_axis.empty()){
            dt = 1.0;
        } else if(_time_axis.size() == 1){
            if(_time_axis[0] > 0.0)
                dt = _time_axis[0];
        } else{
            dt = _time_axis[1] - _time_axis[0];
        }
#if VERBOSE_FIT2X
        std::clog << "-- get_time_resolution: " << dt << std::endl;
#endif
        if(dt == 0.0){
            std::clog << "WARNING dt is zero!" << std::endl;
        }
        return dt;
    }

    int size() const{
        return _size;
    }

    void resize(size_t n);

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
            int scale_start = 0, int scale_stop = -1,
            double scatter_fraction = 0.0,
            double excitation_period = std::numeric_limits<double>::max(),
            double constant_offset = 0.0,
            double number_of_photons=1,
            bool use_amplitude_threshold = false,
            double amplitude_threshold = std::numeric_limits<double>::epsilon(),
            bool use_pile_up_correction = false,
            double instrument_dead_time=120.0,
            double acquisition_time=std::numeric_limits<double>::max(),
            bool add_corrected_irf_as_scatter = false,
            bool scale_model_to_data = false,
            int convolution_method = 0,
            bool take_abs_of_lifetime_spectrum = false,
            bool use_linearization = false,
            double* linearization= nullptr, int n_linearization = -1
    );

    void set_data(double *input, int n_input) {
        set_is_valid(false);
        resize(n_input);
        _data.assign(input, input + n_input);
        set_weights_by_data();
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
            double* data; int n_data;
            get_data(&data, &n_data);
            double re = 0.0;
            for(int i=get_convolution_start(); i<get_convolution_stop(); i++)
                re += data[i];
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
        _convolution_start = v;
    }

    int get_convolution_start() const {
        if(_convolution_start < 0){
            return (int) _irf.size() - _convolution_start;
        } else{
            return _convolution_start;
        }
    }

    void set_convolution_stop(int v) {
        set_is_valid(false);
        _convolution_stop = v;
    }

    int get_convolution_stop() const {
        int stop;
        if(_convolution_stop < 0){
            stop = std::min(size(), std::max(0, (int) size() + _convolution_stop));
        } else if (_convolution_stop == 0){
            stop = size();
        } else{
            stop = _convolution_stop;
        }
        if(stop > size()){
            return size();
        } else{
            return stop;
        }
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
        resize(n_input);
        _irf.assign(input, input + n_input);
        _irf_is_corrected = false;
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
        resize(n_input);
        _linearization_table.assign(input, input + n_input);
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
            for (int i = 0; i < size(); i++)
                _corrected_irf[i] = std::max(_irf[i] - get_irf_background_counts(), 0.0);
            // shift irf
            _corrected_irf = shift_array(
                    _corrected_irf.data(),
                    _corrected_irf.size(),
                    get_irf_shift_channels()
            );
            _irf_is_corrected = true;
        }
        *output_view = _corrected_irf.data();
        *n_output = _corrected_irf.size();
    }

    void get_model(double **output_view, int *n_output){
        if (!get_is_valid()) evaluate();
        *output_view = _model_function.data();
        *n_output = _model_function.size();
    }

    void set_lifetime_spectrum(double *input, int n_input = 0) {
        set_is_valid(false);
        _lifetime_spectrum.assign(input, input+n_input);
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
        resize(n_input);
        _weights.assign(input, input+n_input);
        for(int i=0; i<n_input; i++) _sq_weights[i] = _weights[i] * _weights[i];
    }

    void get_data_weights(double **output_view, int *n_output) {
        *output_view = _weights.data();
        *n_output = _weights.size();
    }

    void set_time_axis(double *input, int n_input) {
        set_is_valid(false);
        resize(n_input);
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
        if(acquisition_time < 0) acquisition_time = std::numeric_limits<double>::max();
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

    void set_tttr_data(
            std::shared_ptr<TTTR> tttr_data, int tttr_micro_time_coarsening
    ){
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
        auto header = tttr_data->get_header();
        double micro_time_resolution = header->get_micro_time_resolution() / tttr_micro_time_coarsening;
        for(size_t i = 0; i < _data.size(); i++) _time_axis.emplace_back(i * micro_time_resolution);
        set_excitation_period(header->get_macro_time_resolution());
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
            double excitation_period = std::numeric_limits<double>::max(),
            bool scale_model_to_data = true, double number_of_photons=1,
            std::vector<double> linearization = std::vector<double>(),
            bool use_linearization = false,
            bool use_pile_up_correction = false, bool use_corrected_irf_as_scatter = false,
            bool use_amplitude_threshold = false,
            double amplitude_threshold = std::numeric_limits<double>::epsilon(),
            double acquisition_time = std::numeric_limits<double>::max(),
            double instrument_dead_time = std::numeric_limits<double>::epsilon(),
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
        for (int i = 0; i < size(); i++)
            _weighted_residuals[i] = (d[i] - m[i]) * w[i];
        *n_output = _weighted_residuals.size();
        *output_view = _weighted_residuals.data();
    }

    static double compute_score(
            double *data, int n_data,
            double *time_axis, int n_time_axis,
            double *irf_histogram, int n_irf_histogram,
            double constant_offset = 0.0,
            double irf_background_counts = 0.0, double scatter_fraction = 0.0,
            double irf_shift_channels = 0.0,
            double *lifetime_spectrum = nullptr, int n_lifetime_spectrum = -1,
            std::vector<int> convolution_range = std::vector<int>({0, -1}),
            std::vector<int> score_range = std::vector<int>({0, -1}),
            double excitation_period = std::numeric_limits<double>::max(),
            bool scale_model_to_data = true, double number_of_photons = 1,
            double *linearization = nullptr, int n_linearization = -1,
            bool use_linearization = false,
            bool use_pile_up_correction = false, bool use_corrected_irf_as_scatter = false,
            bool use_amplitude_threshold = false, double amplitude_threshold = std::numeric_limits<double>::epsilon(),
            double acquisition_time = std::numeric_limits<double>::max(),
            double instrument_dead_time = std::numeric_limits<double>::epsilon(),
            int convolution_method = 0,
            double *data_weights = nullptr, int n_data_weights = -1,
            const char* score_type = "poisson",
            bool take_abs_of_lifetime_spectrum = true
    );

    void evaluate() {
#if VERBOSE_FIT2X
        std::clog << "evaluate..." << std::endl;
#endif
        equalize_length();
        if (!get_is_valid()) {
            double* irf; int n_irf;
            get_corrected_irf(&irf, &n_irf);
            double* scatter; int n_scatter;
            get_irf(&scatter, &n_scatter);
            int s = compute_decay(
                    _model_function.data(), _model_function.size(),
                    _data.data(), _data.size(),
                    _sq_weights.data(),_sq_weights.size(),
                    _time_axis.data(),_time_axis.size(),
                    irf, n_irf,
                    _lifetime_spectrum.data(), _lifetime_spectrum.size(),
                    scatter, n_scatter,
                    get_convolution_start(), get_convolution_stop(),
                    get_score_range()[0], get_score_range()[1],
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
        double mu0 = std::accumulate(decay_histogram.begin(), decay_histogram.end(),0.0);
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
        double dt = get_time_resolution();
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
            const char* score_type= "poisson"
    );

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
    );

};


#endif //FIT2X_DECAY_H

