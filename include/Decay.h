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

#include "ConditionalUpdater.h"
#include "DecayCurve.h"
#include "DecayLifetimeSpectrum.h"
#include "DecayConvolution.h"
#include "DecayScore.h"
#include "DecayBackground.h"
#include "DecayScale.h"
#include "DecayPileup.h"
#include "DecayLinearization.h"

class Decay : public ConditionalUpdater{

public:

    //***********************************************//
    //*      STATIC METHODS                         *//
    //***********************************************//

    /*!
     * Compute a mean lifetime using the moments of the decay and the instrument
     * response function.
     *
     * The lifetime is the first lifetime determined by the method of
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
    ) {
        double m0_irf = std::accumulate(
                irf_histogram.begin(),
                irf_histogram.end(),
                0.0);
        double m1_irf = 0.0;
        for (size_t i = 0; i < irf_histogram.size(); i++)
            m1_irf += i * irf_histogram[i];
        double mu0 = std::accumulate(decay_histogram.begin(), decay_histogram.end(), 0.0);
        double mu1 = 0.0;
        for (size_t i = 0; i < decay_histogram.size(); i++)
            mu1 += i * decay_histogram[i];
        double g1 = mu0 / m0_irf;
        double g2 = (mu1 - g1 * m1_irf) / m0_irf;
        double tau1 = g2 / g1 * micro_time_resolution;
        return tau1;
    }

private:

    std::shared_ptr<DecayScore> decayScore = nullptr;
    std::shared_ptr<DecayLifetimeSpectrum> decayLifetimeSpectrum = nullptr;
    std::shared_ptr<DecayConvolution> decayConvolution = nullptr;
    std::shared_ptr<DecayBackground> decayBackground = nullptr;
    std::shared_ptr<DecayPileup> decayPileup = nullptr;
    std::shared_ptr<DecayScale> decayScale = nullptr;
    std::shared_ptr<DecayLinearization> decayLinearization = nullptr;

    /// The experimental histogram
    std::shared_ptr<DecayCurve> _data = std::make_shared<DecayCurve>();

    /// model function
    std::shared_ptr<DecayCurve> _model = std::make_shared<DecayCurve>();

    /// instrument response function
    std::shared_ptr<DecayCurve> _irf = std::make_shared<DecayCurve>();


public:

    Decay(
            // Data
            std::vector<double> data = std::vector<double>(),
            double acquisition_time = std::numeric_limits<double>::max(),
            std::vector<double> data_weights = std::vector<double>(),
            std::vector<double> time_axis = std::vector<double>(),
            std::vector<double> irf_histogram = std::vector<double>(),
            std::shared_ptr<TTTR> tttr_data = nullptr,
            std::shared_ptr<TTTR> tttr_irf = nullptr,
            int tttr_micro_time_coarsening = 1,
            // Lifetime spectrum
            double* lifetime_spectrum = nullptr, int n_lifetime_spectrum = 0,
            bool use_amplitude_threshold = false,
            bool abs_lifetime_spectrum = false,
            double amplitude_threshold = std::numeric_limits<double>::epsilon(),
            // Convolution
            std::vector<int> convolution_range = std::vector<int>({0, -1}),
            bool use_corrected_irf_as_scatter = true,
            double scatter_fraction = 0.0,
            int convolution_method = ConvFastPeriodicTime,
            double excitation_period = 100., // also for Pile up (repetition rate)
            double irf_shift_channels = 0.0,
            double irf_background_counts = 0,
            // Background
            double constant_offset = 0.0,
            // Pile up
            const char* pile_up_model = "coates",
            double instrument_dead_time = 120,
            bool use_pile_up_correction = false,
            // Decay Scale
            bool scale_model_to_data = false,
            double number_of_photons = -1,
            // Linearization
            double* linearization_table=nullptr, int n_linearization_table=-1,
            bool use_linearization = false,
            // Scoring
            std::vector<int> score_range = std::vector<int>({0, -1}),
            std::string score_type= "poisson"
    );

    ~Decay() = default;

    //***********************************************//
    //*      CONVOLUTION                            *//
    //***********************************************//
    void set_convolution_method(int v){
        set_is_valid(false);
        decayConvolution->set_convolution_method(v);
    }

    int get_convolution_method(){
        return decayConvolution->get_convolution_method();
    }

    void set_excitation_period(double v) {
        set_is_valid(false);
        double rep_rate = 1. / v * 1000.0;
        decayPileup->set_repetition_rate(rep_rate);
        decayConvolution->set_excitation_period(v);
    }

    double get_excitation_period() const {
        return decayConvolution->get_excitation_period();
    }

    void set_irf_shift_channels(double v) {
        set_is_valid(false);
        decayConvolution->set_irf_shift_channels(v);
    }

    double get_irf_shift_channels() const {
        return decayConvolution->get_irf_shift_channels();
    }

    void set_scatter_fraction(double v) {
        set_is_valid(false);
        decayConvolution->set_scatter_fraction(v);
    }

    double get_scatter_fraction() const {
        return decayConvolution->get_scatter_fraction();
    }

    void set_convolution_start(int v) {
        set_is_valid(false);
        decayConvolution->set_convolution_start(v);
        decayScale->set_scale_start(v);
    }

    int get_convolution_start() const {
        return decayConvolution->get_convolution_start();
    }

    void set_convolution_stop(int v) {
        set_is_valid(false);
        decayConvolution->set_convolution_stop(v);
        decayScale->set_scale_stop(v);
    }

    int get_convolution_stop() const {
        return decayConvolution->get_convolution_stop();
    }

    void set_convolution_range(std::vector<int> v) {
        set_convolution_start(v[0]);
        set_convolution_stop(v[1]);
        // set_is_valid(false); // set by set_convolution_start, set_convolution_stop
    }

    std::vector<int> get_convolution_range() const {
        return std::vector<int>({get_convolution_start(), get_convolution_stop()});
    }

    void set_irf(double *input, int n_input) {
#if VERBOSE_FIT2X
        std::clog << "-- set_irf" << std::endl;
#endif
        set_is_valid(false);
        resize(n_input);
        decayConvolution->set_irf(input, n_input);
    }

    void get_irf(double **output_view, int *n_output) {
        decayConvolution->get_irf(output_view, n_output);
    }

    void get_corrected_irf(double **output_view, int *n_output) {
        decayConvolution->get_corrected_irf(output_view, n_output);
    }

    void set_irf_background_counts(double v) {
        set_is_valid(false);
        decayConvolution->set_irf_background_counts(v);
    }

    double get_irf_background_counts() const {
        return decayConvolution->get_irf_background_counts();
    }

    void set_use_corrected_irf_as_scatter(bool v) {
        set_is_valid(false);
        decayConvolution->set_use_corrected_irf_as_scatter(v);
    }

    bool get_use_corrected_irf_as_scatter() const {
        return decayConvolution->get_use_corrected_irf_as_scatter();
    }

    //***********************************************//
    //*      LIFETIME SPECTRUM                      *//
    //***********************************************//
    void set_abs_lifetime_spectrum(bool v){
        set_is_valid(false);
        decayConvolution->set_is_valid(false);
        decayLifetimeSpectrum->set_abs_lifetime_spectrum(v);
    }

    bool get_abs_lifetime_spectrum() const{
        return decayLifetimeSpectrum->get_abs_lifetime_spectrum();
    }

    void set_use_amplitude_threshold(bool v) {
        set_is_valid(false);
        decayConvolution->set_is_valid(false);
        decayLifetimeSpectrum->set_use_amplitude_threshold(v);
    }

    bool get_use_amplitude_threshold() const {
        return decayLifetimeSpectrum->get_use_amplitude_threshold();
    }

    void set_amplitude_threshold(double v) {
        set_is_valid(false);
        decayConvolution->set_is_valid(false);
        decayLifetimeSpectrum->set_amplitude_threshold(v);
    }

    double get_amplitude_threshold() const {
        return decayLifetimeSpectrum->get_amplitude_threshold();
    }

    void set_lifetime_spectrum(double *input, int n_input = 0) {
        set_is_valid(false);
        decayConvolution->set_is_valid(false);
        decayLifetimeSpectrum->set_lifetime_spectrum(input, n_input);
    }

    void get_lifetime_spectrum(double **output_view, int *n_output) {
        decayLifetimeSpectrum->get_lifetime_spectrum(output_view, n_output);
    }

    /*!
     * @brief Get the mean lifetime computed using irf and data
     * 
     * @return double 
     */
    double get_mean_lifetime(){
        return compute_mean_lifetime(
            decayConvolution->corrected_irf->y,
            _data->y,
            micro_time_resolution()
        );
    }

    //***********************************************//
    //*      BACKGROUND                             *//
    //***********************************************//
    void set_constant_offset(double v) {
        set_is_valid(false);
        decayBackground->set_constant_offset(v);
        decayScale->set_constant_background(v);
    }

    double get_constant_offset() const {
        return decayBackground->get_constant_offset();
    }

    //***********************************************//
    //*      PILEUP                                 *//
    //***********************************************//
    void set_use_pile_up_correction(bool v) {
        set_is_valid(false);
        decayPileup->set_use_pile_up_correction(v);
    }

    bool get_use_pile_up_correction() const {
        return decayPileup->get_use_pile_up_correction();
    }

    void set_instrument_dead_time(double v) {
        set_is_valid(false);
        decayPileup->set_instrument_dead_time(v);
    }

    double get_instrument_dead_time() const {
        return decayPileup->get_instrument_dead_time();
    }

    void set_repetition_rate(double v) {
        set_is_valid(false);
        double period = 1. / v * 1000.;
        decayConvolution->set_excitation_period(period);
        decayPileup->set_repetition_rate(v);
    }

    double get_repetition_rate() const {
        return decayPileup->get_repetition_rate();
    }

    void set_pile_up_model(std::string v){
        set_is_valid(false);
        decayPileup->set_pile_up_model(v);
    }

    std::string get_pile_up_model(){
        return decayPileup->get_pile_up_model();
    }

    //***********************************************//
    //*     SCALING                                 *//
    //***********************************************//
    void set_scale_model_to_data(bool v) {
        set_is_valid(false);
        decayScale->set_scale_model_to_data(v);
    }

    bool get_scale_model_to_data() const {
        return decayScale->get_scale_model_to_data();
    }

    void set_number_of_photons(double v) {
        set_is_valid(false);
        decayScale->set_number_of_photons(v);
    }

    double get_number_of_photons(){
        return decayScale->get_number_of_photons();
    }

    //***********************************************//
    //*     LINEARIZATION                           *//
    //***********************************************//
    void set_linearization_table(double *input, int n_input) {
        set_is_valid(false);
        decayLinearization->set_linearization_table(input, n_input);
    }

    void get_linearization_table(double **output_view, int *n_output) {
        decayLinearization->get_linearization_table(output_view, n_output);
    }

    bool get_use_linearization() const{
        return decayLinearization->get_use_linearization();
    }

    void set_use_linearization(bool v){
        set_is_valid(false);
        decayLinearization->set_use_linearization(v);
    }

    /// Computes the model function
    void update_model() {
        decayConvolution->get_decay(*_model);
        decayPileup->add(*_model);
        decayLinearization->add(*_model);
        decayScale->add(*_model);
        decayBackground->add(*_model);
        set_is_valid(true);
    }

    //***********************************************//
    //*     CURVES: EXPERIMENTAL AND MODEL DATA     *//
    //***********************************************//
    double get_time_resolution(){
        return _data->get_dx();
    }

    void set_acquisition_time(double v) {
        _data->set_acquisition_time(v);
    }

    double get_acquisition_time() const {
        return _data->get_acquisition_time();
    }

    void set_data(double *input, int n_input) {
        set_is_valid(false);
        decayConvolution->set_is_valid(false);
        resize(n_input);
        _data->set_y(input, n_input);
    }

    void get_data(double **output_view, int *n_output) {
        *output_view = _data->y.data();
        *n_output = _data->y.size();
    }

    void set_data_weights(double *input, int n_input) {
        set_is_valid(false);
        decayConvolution->set_is_valid(false);
        resize(n_input);
        _data->set_w(input, n_input);
    }

    void get_data_weights(double **output_view, int *n_output) {
        _data->get_w(output_view, n_output);
    }

    int size() const{
        return _data->size();
    }

    void resize(int size){
        _data->resize(size);
        _irf->resize(size);
        _model->resize(size);
        decayLinearization->resize(size);
        decayConvolution->corrected_irf->resize(size);
    }

    void set_time_axis(double *input, int n_input) {
        set_is_valid(false);
        decayConvolution->set_is_valid(false);
        resize(n_input);
        _data->set_x(input, n_input);
        _irf->set_x(input, n_input);
        _model->set_x(input, n_input);
        decayConvolution->corrected_irf->set_x(input, n_input);
    }

    void get_time_axis(double **output_view, int *n_output) {
        _data->get_x(output_view, n_output);
    }

    void get_model(double **output_view, int *n_output){
        if (!get_is_valid())
            update_model();
        _model->get_y(output_view, n_output);
    }

    void set_tttr_data(
            std::shared_ptr<TTTR> tttr,
            int tttr_micro_time_coarsening
    ){
#if VERBOSE_FIT2X
        std::clog << "-- Setting DATA from TTTR..." << std::endl;
#endif
        set_is_valid(false);
        _data->set_tttr(tttr, tttr_micro_time_coarsening);
        set_excitation_period(
                tttr->get_header()->
                get_macro_time_resolution()
                );
    }

    void set_tttr_irf(
            std::shared_ptr<TTTR> tttr,
            int tttr_micro_time_coarsening
    ){
#if VERBOSE_FIT2X
        std::clog << "-- Setting IRF from TTTR..." << std::endl;
#endif
        set_is_valid(false);
        _irf->set_tttr(tttr, tttr_micro_time_coarsening);
        set_excitation_period(
                tttr->get_header()->
                        get_macro_time_resolution()
        );
    }

    double micro_time_resolution(){
        double v = 1.0;
        if(_data->size() > 2){
            v = _data->x[1] - _data->x[0];
        }
        return v;
    }

    //***********************************************//
    //*      SCORE                                  *//
    //***********************************************//
    void set_score_range(std::vector<int> v){
        decayScale->set_scale_start(v[0]);
        decayScale->set_scale_stop(v[1]);
        decayScore->set_score_range(v);
    }

    std::vector<int> get_score_range(){
        return decayScore->get_score_range();
    }

    void get_weighted_residuals(double **output_view, int *n_output) {
        if(!get_is_valid()){
            update_model();
        }
        decayScore->get_weighted_residuals(output_view, n_output);
    }

    double get_score(
            int x_min = -1,
            int x_max = -1,
            const char* score_type= "poisson"
    ){
        if(!get_is_valid()){
            update_model();
        }
        return decayScore->get_score(x_min, x_max, score_type);
    }

    std::string get_score_type(){
        return decayScore->get_score_type();
    }

    void set_score_type(std::string v){
        decayScore->set_score_type(v);
    }


};


#endif //FIT2X_DECAY_H

