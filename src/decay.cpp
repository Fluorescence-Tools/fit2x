#include "decay.h"

std::vector<double> Decay::shift_array(
        double* input, int n_input,
        double shift,
        bool set_outside,
        double outside_value
){
    auto out = std::vector<double>(n_input);
    if (n_input > 0) {
        // mod
        int r = (int) shift % n_input;
        int ts_i = r < 0 ? r + n_input : r;
        double ts_f = shift - std::floor(shift);
#if VERBOSE_FIT2X
        std::clog << "shift_array..." << std::endl;
        std::clog << "-- n_input: " << n_input << std::endl;
        std::clog << "-- shift: " << shift << std::endl;
        std::clog << "-- shift integer: " << ts_i << std::endl;
        std::clog << "-- shift float: " << ts_f << std::endl;
#endif
        std::vector<double> tmp1(input, input + n_input);
        std::vector<double> tmp2(input, input + n_input);
        std::rotate(tmp1.begin(), tmp1.begin() + ts_i, tmp1.end());
        std::rotate(tmp2.begin(), tmp2.begin() + (ts_i + 1), tmp2.end());
        for (int i = 0; i < n_input; i++) {
            out[i] = tmp1[i] * (1.0 - ts_f) + tmp2[i] * ts_f;
        }
        if (set_outside) {
            if (shift > 0) {
                for (int i = 0; i < r; i++)
                    out[i] = outside_value;
            } else if (shift < 0) {
                for (int i = n_input - r; i < n_input; i++)
                    out[i] = outside_value;
            }
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
    *n_output = std::min(n_curve1, n_curve2);
    start = std::min(0, start);
    stop = stop < 0? *n_output: std::min(*n_output, stop);
    double sum_curve_1 = std::accumulate( curve1 + start, curve1 + stop, 0.0);
    double sum_curve_2 = std::accumulate( curve2 + start, curve2 + stop, 0.0);
    double f1 = (1. - areal_fraction_curve2);
    double f2 = areal_fraction_curve2 * sum_curve_1 / sum_curve_2;
#ifndef _WIN32
#pragma omp simd
#endif
    for(int i=start; i<stop;i++) curve1[i] *= f1;
#ifndef _WIN32
#pragma omp simd
#endif
    for(int i=start; i<stop;i++) curve1[i] += f2 * curve2[i];
    *output = curve1;
    *n_output = n_curve1;
}


void Decay::resize(size_t n){
#if VERBOSE_FIT2X
    std::clog << "-- Decay::resize - new size: " << n << std::endl;
#endif
    _irf.resize(n, 0.0);
    _corrected_irf.resize(n, 0.0);
    _model_function.resize(n, 0.0);
    _data.resize(n, 0.0);
    _weights.resize(n, 0.0);
    _sq_weights.resize(n, 0.0);
    _linearization_table.resize(n, 1.0);
    _weighted_residuals.resize(n, 0.0);
    if(_time_axis.size() < n){
        double dt = get_time_resolution();
        size_t start = _time_axis.size();
        _time_axis.resize(n);
        for(size_t i=start; i < _time_axis.size(); i++)
            _time_axis[i] = _time_axis[i - 1] + dt;
    }
    _size = n;
}


double Decay::get_score(int x_min, int x_max, const char* score_type){
#if VERBOSE_FIT2X
    std::clog << "CHI2" << std::endl;
    std::clog << "-- data range: " << x_min << ", " << x_max << std::endl;
#endif
    x_min = (x_min < 0) ? std::max(this->_score_range_min, 0) : x_min;
    x_max = (x_max < 0) ? (int)_data.size(): std::min(x_max, (int)_data.size());
    double v = 0.0;
    if(strcmp(score_type, "normal") == 0) {
        double* wres; int nwres;
        get_weighted_residuals(&wres, &nwres);
        auto sr = get_score_range();
        for(int i=sr[0];i<sr[1];i++) v += wres[i] * wres[i];
    } else{
        double* m; int nm; get_model(&m, &nm);
        double* d; int nd; get_data(&d, &nd);
        auto data = std::vector<double>(d, d+nd);
        auto model = std::vector<double>(m, m+nm);
        v = statistics::chi2_counting(data, model, x_min, x_max, score_type);
    }
#if VERBOSE_FIT2X
    std::clog << "-- x_min: " << x_min << std::endl;
        std::clog << "-- x_max: " << x_max << std::endl;
        std::clog << "-- chi2: " << v << std::endl;
#endif
    if(get_is_valid()) return v;
    else return INFINITY;
}

double Decay::compute_score(
        double *data, int n_data,
        double *time_axis, int n_time_axis,
        double *irf_histogram, int n_irf_histogram,
        double constant_offset,
        double irf_background_counts, double scatter_fraction,
        double irf_shift_channels,
        double *lifetime_spectrum, int n_lifetime_spectrum,
        std::vector<int> convolution_range,
        std::vector<int> score_range,
        double excitation_period,
        bool scale_model_to_data, double number_of_photons,
        double *linearization, int n_linearization,
        bool use_linearization,
        bool use_pile_up_correction, bool use_corrected_irf_as_scatter,
        bool use_amplitude_threshold, double amplitude_threshold,
        double acquisition_time, double instrument_dead_time,
        int convolution_method,
        double *data_weights, int n_data_weights,
        const char* score_type,
        bool take_abs_of_lifetime_spectrum
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
            model_function.data(), model_function.size(),
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
#if VERBOSE_FIT2X
        std::clog << "COMPUTING SCORE..." << std::endl;
            std::clog << "-- score_range: " << score_range[0] << std::endl;
            std::clog << "-- score_range: " << score_range[1] << std::endl;
            std::clog << "-- score_type: " << score_type << std::endl;
#endif
        double score = statistics::chi2_counting(
                decay_histogram,
                model_function,
                score_range[0], score_range[1],
                score_type
        );
        return score;
    } else{
        return std::numeric_limits<double>::infinity();
    }
}


void Decay::set(
        std::vector<double> data,
        std::vector<double> time_axis,
        std::vector<double> irf_histogram,
        double constant_offset,
        double irf_background_counts,
        double scatter_fraction,
        double irf_shift_channels,
        std::vector<double> lifetime_spectrum,
        std::vector<int> convolution_range,
        std::vector<int> score_range,
        double excitation_period,
        bool scale_model_to_data,
        double number_of_photons,
        std::vector<double> linearization,
        bool use_linearization,
        bool use_pile_up_correction,
        bool use_corrected_irf_as_scatter,
        bool use_amplitude_threshold,
        double amplitude_threshold,
        double acquisition_time,
        double instrument_dead_time,
        int convolution_method,
        std::vector<double> data_weights,
        std::shared_ptr<TTTR> tttr_data,
        std::shared_ptr<TTTR> tttr_irf,
        int tttr_micro_time_coarsening,
        bool force_fill_irf,
        bool force_time_axis
) {

    // set data
    if (tttr_data != nullptr) {
        set_tttr_data(tttr_data, tttr_micro_time_coarsening);
        // sets also excitation period
    } else{
        if(!data.empty()){
            set_data(data.data(), data.size());
            if(data_weights.empty()) set_weights_by_data();
        }
    }

    // set time axis
    if(!time_axis.empty()){
        set_time_axis(time_axis.data(), time_axis.size());
    } else if(force_time_axis){
        auto t = std::vector<double>(size());
        for(int i=0; i < size(); i++) t[i] = i;
        set_time_axis(t.data(), t.size());
    }

    // set irf
    if (tttr_irf != nullptr)
        set_tttr_irf(tttr_irf, tttr_micro_time_coarsening);
    else if(!irf_histogram.empty()){
        set_irf(irf_histogram.data(), irf_histogram.size());
    }
    else if(force_fill_irf){
        if(_irf.empty()) resize(1);
        _irf[0] = 1.0;
    }

    // set linearization
    if(!linearization.empty()){
        set_linearization(linearization.data(), linearization.size());
        set_use_linearization(use_linearization);
    } else{
        set_use_linearization(false);
    }

    // set weights
    if(!data_weights.empty()){
        set_data_weights(data_weights.data(), data_weights.size());
    }

    // make sure that score range is of correct length
    if (score_range.empty()) {
        score_range.emplace_back(0);
        score_range.emplace_back(-1);
    } else if (score_range.size() < 2) {
        score_range.emplace_back(-1);
    }

    // convolution method
    set_convolution_method(convolution_method);
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
}


void Decay::scale_model(
        bool scale_model_to_data,
        double number_of_photons,
        int start, int stop,
        double constant_background,
        double* model,
        double* data,
        double* squared_data_weights
        ){
    double scale = 1.0;
    if(scale_model_to_data || (number_of_photons < 0)){
        // scale the area to the data in the range start, stop
#if VERBOSE_FIT2X
        std::clog << "-- scaling model to data..." << std::endl;
#endif
        rescale_w_bg(
                model, data, squared_data_weights, constant_background,
                &scale, start, stop
        );
    } else{
#if VERBOSE_FIT2X
        std::clog << "-- scaling decay to " << number_of_photons << " photons." << std::endl;
#endif
        // normalize model to total_area
        double number_of_photons_model = std::accumulate(
                model + start,
                model + stop,
                0.0
        );
        scale = number_of_photons / number_of_photons_model;
        for(int i=start; i < stop; i++) model[i] *= scale;
    }
}

int Decay::compute_decay(
        double* model_function, int n_model_function,
        double* data, int n_data,
        double* squared_data_weights, int n_squared_data_weights,
        double* time_axis, int n_time_axis,
        double* irf_histogram, int n_irf_histogram,
        double* lifetime_spectrum, int n_lifetime_spectrum,
        double* scatter, int n_scatter,
        int convolution_start, int convolution_stop,
        double scatter_fraction,
        double excitation_period,
        double constant_offset,
        double number_of_photons,
        bool use_amplitude_threshold,
        double amplitude_threshold,
        bool use_pile_up_correction,
        double instrument_dead_time,
        double acquisition_time,
        bool use_corrected_irf_as_scatter,
        bool scale_model_to_data,
        int convolution_method,
        bool take_abs_of_lifetime_spectrum,
        bool use_linearization,
        double* linearization, int n_linearization
){
#if VERBOSE_FIT2X
    std::clog << "compute_decay..." << std::endl;
    std::clog << "-- n_model_function: " << n_model_function << std::endl;
    std::clog << "-- n_data: " << n_data << std::endl;
    std::clog << "-- n_time_axis: " << n_time_axis << std::endl;
    std::clog << "-- n_instrument_response_function: " << n_irf_histogram << std::endl;
    std::clog << "-- n_lifetime_spectrum: " << n_lifetime_spectrum << std::endl;
    std::clog << "-- convolution_start: " << convolution_start << std::endl;
    std::clog << "-- convolution_stop: " << convolution_stop << std::endl;
    std::clog << "-- scatter_fraction: " << scatter_fraction << std::endl;
    std::clog << "-- excitation_period: " << excitation_period << std::endl;
    std::clog << "-- constant_offset: " << constant_offset << std::endl;
    std::clog << "-- number_of_photons: " << number_of_photons << std::endl;
    std::clog << "-- use_amplitude_threshold: " << use_amplitude_threshold << std::endl;
    std::clog << "-- amplitude_threshold: " << amplitude_threshold << std::endl;
    std::clog << "-- use_pile_up_correction: " << use_pile_up_correction << std::endl;
    std::clog << "-- instrument_dead_time: " << instrument_dead_time << std::endl;
    std::clog << "-- acquisition_time: " << acquisition_time << std::endl;
    std::clog << "-- use_corrected_irf_as_scatter: " << use_corrected_irf_as_scatter << std::endl;
    std::clog << "-- use_pile_up_correction: " << use_pile_up_correction << std::endl;
    std::clog << "-- use_pile_up_correction: " << use_pile_up_correction << std::endl;
    std::clog << "-- scale_model_to_data: " << scale_model_to_data << std::endl;
    std::clog << "-- convolution_method: " << convolution_method << std::endl;
    std::clog << "-- take_abs_of_lifetime_spectrum: " << take_abs_of_lifetime_spectrum << std::endl;
    std::clog << "-- use_linearization: " << use_linearization << std::endl;
    std::clog << "-- use_corrected_irf_as_scatter: " << use_corrected_irf_as_scatter << std::endl;
#endif
    if(use_amplitude_threshold){
        discriminate_small_amplitudes(
                lifetime_spectrum, n_lifetime_spectrum,
                amplitude_threshold
        );
    }
    // Take abs(lifetime_spectrum)
    std::vector<double> lt; lt.resize(n_lifetime_spectrum);
    if(take_abs_of_lifetime_spectrum) {
        for(int i=0; i<n_lifetime_spectrum; i++) lt[i] = std::abs(lifetime_spectrum[i]);
    } else{
        for(int i=0; i<n_lifetime_spectrum; i++) lt[i] = lifetime_spectrum[i];
    }
    if(convolution_method == 0){
        fconv_per_cs_time_axis(
                model_function, n_model_function,
                time_axis, n_time_axis,
                irf_histogram, n_irf_histogram,
                lt.data(), lt.size(),
                convolution_start, convolution_stop,
                excitation_period
        );
    } else if(convolution_method == 1){
        fconv_cs_time_axis(
                model_function, n_model_function,
                time_axis, n_time_axis,
                irf_histogram, n_irf_histogram,
                lt.data(), lt.size(),
                convolution_start, convolution_stop
        );
    } else if (convolution_method == 2) {
        double dt = time_axis[1] - time_axis[0];
        fconv_per(
                model_function, lt.data(),
                irf_histogram, (int) lt.size() / 2,
                convolution_start, convolution_stop,
                n_model_function, excitation_period,
                dt
        );
    } else if (convolution_method == 3) {
        double dt = time_axis[1] - time_axis[0];
        fconv(
                model_function, lt.data(),
                irf_histogram, (int) lt.size() / 2,
                convolution_start, n_model_function,
                dt
        );
    } else if (convolution_method == 4) {
        double dt = time_axis[1] - time_axis[0];
        fconv_avx(
                model_function, lt.data(),
                irf_histogram, (int) lt.size() / 2,
                convolution_start, n_model_function,
                dt
        );
    }  else if (convolution_method == 5) {
        double dt = time_axis[1] - time_axis[0];
        fconv_per_avx(
                model_function, lt.data(),
                irf_histogram, (int) lt.size() / 2,
                convolution_start, convolution_stop,
                n_model_function, excitation_period,
                dt
        );
    }
#if VERBOSE_FIT2X
    std::clog << "model_function [:64]: ";
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
#if VERBOSE_FIT2X
    std::clog << "model_with_scatter [:64]: ";
    for(int i=0; i<64; i++) std::clog << model_function[i] << " ";
    std::clog << std::endl;
#endif
    if(use_pile_up_correction){
        double rep_rate = 1. / excitation_period * 1000.;
        add_pile_up_to_model(
                model_function, n_model_function,
                data, n_data,
                rep_rate, instrument_dead_time,
                acquisition_time
        );
    }
#if VERBOSE_FIT2X
    std::clog << "model_with_scatter - pile up [:64]: ";
    for(int i=0; i<64; i++) std::clog << model_function[i] << " ";
    std::clog << std::endl;
#endif
    // linearization
    if(use_linearization && (n_linearization > 0) && (linearization!= nullptr)){
        for(int i = 0; i < std::min(n_linearization, n_model_function); i++){
            model_function[i] *= linearization[i];
        }
    }
    // scale model function
    scale_model(
            scale_model_to_data,
            number_of_photons,
            convolution_start, convolution_stop,
            constant_offset, model_function, data, squared_data_weights
    );
#if VERBOSE_FIT2X
    std::clog << "model_with_scatter - scaled [:64]: ";
    for(int i=0; i<64; i++) std::clog << model_function[i] << " ";
    std::clog << std::endl;
#endif
#if VERBOSE_FIT2X
    std::clog << "Adding Background" << std::endl;
    std::clog << "-- add constant offset: " << constant_offset << std::endl;
#endif
    for(int i=0; i<n_model_function;i++)
        model_function[i] += constant_offset;
#if VERBOSE_FIT2X
    std::clog << "model final [:64]: ";
    for(int i=0; i<64; i++) std::clog << model_function[i] << " ";
    std::clog << std::endl;
#endif
    return 0;
}

