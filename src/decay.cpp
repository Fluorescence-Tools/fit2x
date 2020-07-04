#include "decay.h"



void Decay::shift_array(
        double* input, int n_input,
        double** output, int *n_output,
        double shift,
        bool set_outside,
        double outside_value
){
    // mod
    int r = (int) shift % n_input;
    int ts_i = r < 0 ? r + n_input : r;
    double ts_f = shift - std::floor(shift);
#if VERBOSE
    std::clog << "shift_array..." << std::endl;
    std::clog << "-- n_input: " << n_input << std::endl;
    std::clog << "-- shift: " << shift << std::endl;
    std::clog << "-- shift integer: " << ts_i << std::endl;
    std::clog << "-- shift float: " << ts_f << std::endl;
#endif
    auto tmp = (double*) calloc(n_input, sizeof(double));
    std::vector<double> tmp1(input, input+n_input);
    std::vector<double> tmp2(input, input+n_input);
    std::rotate(tmp1.begin(), tmp1.begin()+ts_i, tmp1.end());
    std::rotate(tmp2.begin(), tmp2.begin()+(ts_i + 1), tmp2.end());
    for(int i=0; i < n_input; i++){
        tmp[i] = tmp1[i] * (1.0 - ts_f) + tmp2[i] * ts_f;
    }
    if(set_outside){
        if(shift > 0){
            for(int i=0; i<r; i++)
                tmp[i] = outside_value;
        } else if(shift<0){
            for(int i = n_input - r; i<n_input; i++)
                tmp[i] = outside_value;
        }
    }
    *output = tmp;
    *n_output = n_input;
}


void Decay::add_curve(
        double** output, int *n_output,
        double* curve1, int n_curve1,
        double* curve2, int n_curve2,
        double areal_fraction_curve2,
        int start,
        int stop
){
#if VERBOSE
    std::clog << "add_curve..." << std::endl;
    std::clog << "-- start: " << start << std::endl;
    std::clog << "-- stop: " << stop << std::endl;
    std::clog << "-- areal_fraction_curve2: " << areal_fraction_curve2 << std::endl;
#endif
    *n_output = std::min(n_curve1, n_curve2);
    start = std::min(0, start);
    stop = stop < 0? *n_output: std::min(*n_output, stop);
    auto tmp  = (double*) malloc(*n_output * sizeof(double));
    auto tmp1 = (double*) malloc(n_curve1 * sizeof(double));
    auto tmp2 = (double*) malloc(n_curve2 * sizeof(double));
    double sum_curve_1 = std::accumulate( curve1 + start, curve1 + stop, 0.0);
    double sum_curve_2 = std::accumulate( curve2 + start, curve2 + stop, 0.0);
    for(int i=0; i<n_curve1;i++) tmp1[i]  = curve1[i] / sum_curve_1;
    for(int i=0; i<n_curve2;i++) tmp2[i]  = curve2[i] / sum_curve_2;
    for(int i=start; i<stop;i++)
        tmp[i] = (
                    (1. - areal_fraction_curve2) * tmp1[i] + areal_fraction_curve2 * tmp2[i]
                ) * sum_curve_1;
    *output = tmp;
}


void Decay::scale_model(
        bool scale_model_to_data,
        double number_of_photons,
        int convolution_start, int convolution_stop,
        double constant_background,
        double* model,
        double* data,
        double* squared_data_weights
        ){
    double scale = 0.0;
    if(scale_model_to_data || (number_of_photons < 0)){
        // scale the area to the data in the range start, stop
#if VERBOSE
        std::clog << "-- scaling model to data..." << std::endl;
#endif
        rescale_w_bg(
                model,
                data, squared_data_weights,
                constant_background,
                &scale, convolution_start, convolution_stop
        );
    } else{
#if VERBOSE
        std::clog << "-- scaling decay to " << number_of_photons << " photons." << std::endl;
#endif
        // normalize model to total_area
        double number_of_photons_model = std::accumulate(
                model + convolution_start,
                model + convolution_stop,
                0.0
        );
        scale = number_of_photons / number_of_photons_model;
        for(int i=convolution_start; i < convolution_stop; i++) model[i] *= scale;
    }
}

void Decay::compute_decay(
        double* model_function, int n_model_function,
        double* data, int n_data,
        double* squared_weights, int n_weights,
        double* time_axis, int n_time_axis,
        double* instrument_response_function, int n_instrument_response_function,
        double* lifetime_spectrum, int n_lifetime_spectrum,
        int convolution_start, int convolution_stop,
        double irf_background_counts,
        double irf_shift_channels,
        double scatter_fraction,
        double excitation_period,
        double constant_background,
        double number_of_photons,
        bool use_amplitude_threshold,
        double amplitude_threshold,
        bool add_pile_up,
        double instrument_dead_time,
        double acquisition_time,
        bool use_corrected_irf_as_scatter,
        bool scale_model_to_data
){
    convolution_stop = convolution_stop > 0 ?
                       std::min({n_time_axis, n_instrument_response_function, n_model_function, convolution_stop}) :
                       std::min({n_time_axis, n_instrument_response_function, n_model_function});
    convolution_stop -= 1;
#if VERBOSE
    std::clog << "compute_decay..." << std::endl;
    std::clog << "-- n_model_function: " << n_model_function << std::endl;
    std::clog << "-- n_data: " << n_data << std::endl;
    std::clog << "-- n_weights: " << n_weights << std::endl;
    std::clog << "-- n_time_axis: " << n_time_axis << std::endl;
    std::clog << "-- n_instrument_response_function: " << n_instrument_response_function << std::endl;
    std::clog << "-- n_lifetime_spectrum: " << n_lifetime_spectrum << std::endl;
    std::clog << "-- convolution_start: " << convolution_start << std::endl;
    std::clog << "-- convolution_stop: " << convolution_stop << std::endl;
    std::clog << "-- scale_model_to_data: " << scale_model_to_data << std::endl;
    std::clog << "-- irf_background_counts: " << irf_background_counts << std::endl;
    std::clog << "-- irf_shift_channels: " << irf_shift_channels << std::endl;
    std::clog << "-- irf_areal_fraction: " << scatter_fraction << std::endl;
    std::clog << "-- period: " << excitation_period << std::endl;
    std::clog << "-- constant_background: " << constant_background << std::endl;
    std::clog << "-- number_of_photons: " << number_of_photons << std::endl;
    std::clog << "-- use_amplitude_threshold: " << use_amplitude_threshold << std::endl;
    std::clog << "-- amplitude_threshold: " << amplitude_threshold << std::endl;
    std::clog << "-- add_pile_up: " << add_pile_up << std::endl;
    std::clog << "-- use_corrected_irf_as_scatter: " << use_corrected_irf_as_scatter << std::endl;
#endif
    // correct irf for background counts
    auto irf_bg_corrected = std::vector<double>(n_instrument_response_function);
    for(int i=0; i < n_instrument_response_function; i++)
        irf_bg_corrected[i] = std::max(
                0.0, instrument_response_function[i] - irf_background_counts
        );

    // shift irf
    double* irf_bg_shift_corrected; int n_irf_bg_shift_corrected;
    shift_array(
        irf_bg_corrected.data(), irf_bg_corrected.size(),
        &irf_bg_shift_corrected, &n_irf_bg_shift_corrected,
        irf_shift_channels
    );
    // convolve lifetime spectrum with irf
    fconv_per_cs_time_axis(
            model_function, n_model_function,
            time_axis, n_time_axis,
            irf_bg_shift_corrected, n_irf_bg_shift_corrected,
            lifetime_spectrum, n_lifetime_spectrum,
            convolution_start, convolution_stop,
            use_amplitude_threshold, amplitude_threshold,
            excitation_period
    );
    free(irf_bg_shift_corrected);

    // add scatter fraction (irf)
    double* model_with_scatter; int n_decay_irf;
    if(use_corrected_irf_as_scatter){
        add_curve(
                &model_with_scatter, &n_decay_irf,
                model_function, n_model_function,
                irf_bg_shift_corrected, n_irf_bg_shift_corrected,
                scatter_fraction,
                convolution_start, convolution_stop
        );
    } else{
        add_curve(
                &model_with_scatter, &n_decay_irf,
                model_function, n_model_function,
                instrument_response_function, n_instrument_response_function,
                scatter_fraction,
                convolution_start, convolution_stop
        );
    }

    if(add_pile_up){
        double rep_rate = 1. / excitation_period * 1000.;
        add_pile_up_to_model(
                model_with_scatter, n_decay_irf,
                data, n_data,
                rep_rate, instrument_dead_time,
                acquisition_time
        );
    }

    // scale model function
    scale_model(
            scale_model_to_data,
            number_of_photons,
            convolution_start, convolution_stop,
            constant_background, model_with_scatter, data, squared_weights
    );

#if VERBOSE
    std::clog << "Adding Background" << std::endl;
    std::clog << "-- add constant background: " << constant_background << std::endl;
#endif
    for(int i=0; i<n_model_function;i++)
        model_function[i] = model_with_scatter[i] + constant_background;
}

