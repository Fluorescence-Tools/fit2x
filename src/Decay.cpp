#include "Decay.h"


Decay::Decay(
        // Data
        std::vector<double> data,
        double acquisition_time,
        std::vector<double> data_weights,
        std::vector<double> time_axis,
        std::vector<double> irf_histogram,
        std::shared_ptr<TTTR> tttr_data,
        std::shared_ptr<TTTR> tttr_irf,
        int tttr_micro_time_coarsening,
        // Lifetime spectrum
        double *lifetime_spectrum, int n_lifetime_spectrum,
        bool use_amplitude_threshold,
        bool abs_lifetime_spectrum,
        double amplitude_threshold,
        // Convolution
        std::vector<int> convolution_range,
        bool use_corrected_irf_as_scatter,
        double scatter_fraction,
        int convolution_method,
        double excitation_period,
        double irf_shift_channels,
        double irf_background_counts,
        // Background
        double constant_offset,
        // Pile up
        const char *pile_up_model,
        double instrument_dead_time,
        bool use_pile_up_correction,
        // Decay Scale
        bool scale_model_to_data,
        double number_of_photons,
        // Linearization
        std::vector<double> linearization_table,
        bool use_linearization,
        // Scoring
        std::vector<int> score_range,
        std::string score_type
) {
    decayScore = new DecayScore(_model, &_data, score_range, score_type);
    decayLifetimeSpectrum = new DecayLifetimeSpectrum(
            lifetime_spectrum, n_lifetime_spectrum,
            use_amplitude_threshold,
            abs_lifetime_spectrum,
            amplitude_threshold);
    decayConvolution = new DecayConvolution(
            &_irf,
            decayLifetimeSpectrum,
            convolution_range,
            use_corrected_irf_as_scatter,
            scatter_fraction,
            convolution_method,
            excitation_period,
            irf_shift_channels,
            irf_background_counts);
    _model = &decayConvolution->decay;
    decayBackground = new DecayBackground(constant_offset);
    double repetition_rate = 1. / excitation_period *  1000.0;
    decayPileup = new DecayPileup(&_data, pile_up_model,
                                  repetition_rate,
                                  instrument_dead_time,
                                  use_pile_up_correction);
    decayScale = new DecayScale(&_data,
                                scale_model_to_data,
                                number_of_photons,
                                score_range,
                                constant_offset
    );
    decayLinearization = new DecayLinearization();
    set_linearization_table(linearization_table.data(), linearization_table.size());
    set_use_linearization(use_linearization);
    decayScore = new DecayScore(_model, &_data, score_range, score_type);
//
//    _data.set_tttr(tttr_data, tttr_micro_time_coarsening);
//    _irf.set_tttr(tttr_irf, tttr_micro_time_coarsening);
//    if(!time_axis.empty())
//        set_time_axis(time_axis.data(), time_axis.size());
//    if(!irf_histogram.empty())
//        set_irf(irf_histogram.data(), irf_histogram.size());
//    if(!data.empty())
//        set_data(data.data(), data.size());
//    if(!data_weights.empty())
//        set_data_weights(data_weights.data(), data_weights.size());
}
