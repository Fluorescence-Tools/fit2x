#include "decay.h"



std::vector<double> Decay::shift_array(
        double* input, int n_input,
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
    auto out = std::vector<double>(n_input);
    std::vector<double> tmp1(input, input+n_input);
    std::vector<double> tmp2(input, input+n_input);
    std::rotate(tmp1.begin(), tmp1.begin()+ts_i, tmp1.end());
    std::rotate(tmp2.begin(), tmp2.begin()+(ts_i + 1), tmp2.end());
    for(int i=0; i < n_input; i++){
        out[i] = tmp1[i] * (1.0 - ts_f) + tmp2[i] * ts_f;
    }
    if(set_outside){
        if(shift > 0){
            for(int i=0; i<r; i++)
                out[i] = outside_value;
        } else if(shift<0){
            for(int i = n_input - r; i<n_input; i++)
                out[i] = outside_value;
        }
    }
    return out;
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
    double sum_curve_1 = std::accumulate( curve1 + start, curve1 + stop, 0.0);
    double sum_curve_2 = std::accumulate( curve2 + start, curve2 + stop, 0.0);
    double f1 = (1. - areal_fraction_curve2);
    double f2 = areal_fraction_curve2 * sum_curve_1 / sum_curve_2;
    for(int i=start; i<stop;i++) curve1[i] = (f1 * curve1[i] + f2 * curve2[i]);
    *output = curve1;
    *n_output = n_curve1;
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
        double* irf_histogram, int n_irf_histogram,
        double* lifetime_spectrum, int n_lifetime_spectrum,
        double* scatter, int n_scatter,
        int convolution_start, int convolution_stop,
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
        bool scale_model_to_data,
        int convolution_method
){
    convolution_stop = convolution_stop > 0 ?
                       std::min({n_time_axis, n_irf_histogram, n_model_function, convolution_stop}) :
                       std::min({n_time_axis, n_irf_histogram, n_model_function});
    // decremente to make sure to stay in bounds
    convolution_stop = convolution_stop - 1;
#if VERBOSE
    std::clog << "compute_decay..." << std::endl;
    std::clog << "-- n_model_function: " << n_model_function << std::endl;
    std::clog << "-- n_data: " << n_data << std::endl;
    std::clog << "-- n_weights: " << n_weights << std::endl;
    std::clog << "-- n_time_axis: " << n_time_axis << std::endl;
    std::clog << "-- n_instrument_response_function: " << n_irf_histogram << std::endl;
    std::clog << "-- n_lifetime_spectrum: " << n_lifetime_spectrum << std::endl;
    std::clog << "-- convolution_start: " << convolution_start << std::endl;
    std::clog << "-- convolution_stop: " << convolution_stop << std::endl;
    std::clog << "-- scale_model_to_data: " << scale_model_to_data << std::endl;
    std::clog << "-- irf_areal_fraction: " << scatter_fraction << std::endl;
    std::clog << "-- period: " << excitation_period << std::endl;
    std::clog << "-- constant_background: " << constant_background << std::endl;
    std::clog << "-- number_of_photons: " << number_of_photons << std::endl;
    std::clog << "-- use_amplitude_threshold: " << use_amplitude_threshold << std::endl;
    std::clog << "-- amplitude_threshold: " << amplitude_threshold << std::endl;
    std::clog << "-- add_pile_up: " << add_pile_up << std::endl;
    std::clog << "-- use_corrected_irf_as_scatter: " << use_corrected_irf_as_scatter << std::endl;
#endif
    if(use_amplitude_threshold){
        discriminate_small_amplitudes(
                lifetime_spectrum, n_lifetime_spectrum,
                amplitude_threshold
        );
    }
    if(convolution_method == 0){
        fconv_per_cs_time_axis(
                model_function, n_model_function,
                time_axis, n_time_axis,
                irf_histogram, n_irf_histogram,
                lifetime_spectrum, n_lifetime_spectrum,
                convolution_start, convolution_stop,
                excitation_period
        );
    } else if(convolution_method == 1){
        fconv_cs_time_axis(
                model_function, n_model_function,
                time_axis, n_time_axis,
                irf_histogram, n_irf_histogram,
                lifetime_spectrum, n_lifetime_spectrum,
                convolution_start, convolution_stop
        );
    } else if (convolution_method == 2) {
        double dt = time_axis[1] - time_axis[0];
        fconv_per(
                model_function, lifetime_spectrum,
                irf_histogram, n_lifetime_spectrum,
                convolution_start, convolution_stop,
                n_model_function, excitation_period,
                dt
        );
    } else if (convolution_method == 3) {
        double dt = time_axis[1] - time_axis[0];
        fconv(
                model_function, lifetime_spectrum,
                irf_histogram, n_lifetime_spectrum,
                convolution_start, convolution_stop,
                dt
        );
    }
#if VERBOSE
    std::clog << "fconv_per_cs_time_axis - model_function [:64]: ";
    for(int i=0; i<64; i++) std::clog << model_function[i] << " ";
    std::clog << std::endl;
#endif
    // add scatter fraction (irf)
    if(use_corrected_irf_as_scatter || (scatter == nullptr) || (n_scatter <= 0)){
        add_curve(
                &model_function, &n_model_function,
                model_function, n_model_function,
                irf_histogram, n_irf_histogram,
                scatter_fraction,
                convolution_start, convolution_stop
        );
    } else{
        add_curve(
                &model_function, &n_model_function,
                model_function, n_model_function,
                scatter, n_scatter,
                scatter_fraction,
                convolution_start, convolution_stop
        );
    }
#if VERBOSE
    std::clog << "model_with_scatter [:64]: ";
    for(int i=0; i<64; i++) std::clog << model_function[i] << " ";
    std::clog << std::endl;
#endif
    if(add_pile_up){
        double rep_rate = 1. / excitation_period * 1000.;
        add_pile_up_to_model(
                model_function, n_model_function,
                data, n_data,
                rep_rate, instrument_dead_time,
                acquisition_time
        );
    }
#if VERBOSE
    std::clog << "model_with_scatter - pile up [:64]: ";
    for(int i=0; i<64; i++) std::clog << model_function[i] << " ";
    std::clog << std::endl;
#endif
    // scale model function
    scale_model(
            scale_model_to_data,
            number_of_photons,
            convolution_start, convolution_stop,
            constant_background, model_function, data, squared_weights
    );
#if VERBOSE
    std::clog << "model_with_scatter - scaled [:64]: ";
    for(int i=0; i<64; i++) std::clog << model_function[i] << " ";
    std::clog << std::endl;
#endif
#if VERBOSE
    std::clog << "Adding Background" << std::endl;
    std::clog << "-- add constant background: " << constant_background << std::endl;
#endif
    for(int i=0; i<n_model_function;i++)
        model_function[i] += constant_background;
#if VERBOSE
    std::clog << "model final [:64]: ";
    for(int i=0; i<64; i++) std::clog << model_function[i] << " ";
    std::clog << std::endl;
#endif
}

