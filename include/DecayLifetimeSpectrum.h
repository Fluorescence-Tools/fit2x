#ifndef FIT2X_DECAYLIFETIMESPECTRUM_H
#define FIT2X_DECAYLIFETIMESPECTRUM_H

#include <vector>
#include <limits>

#include "fsconv.h"

class DecayLifetimeSpectrum{

private:

    /// Lifetime spectrum / original
    std::vector<double> _lifetime_spectrum = std::vector<double>();

    /// Lifetime spectrum / for getter
    std::vector<double> lt = std::vector<double>();

    /// Threshold used to discriminate lifetimes with small amplitudes
    double amplitude_threshold = std::numeric_limits<double>::epsilon();

    /// If true lifetimes with small amplitudes are discriminated
    bool use_amplitude_threshold = false;

    /// If true absolute values of lifetime spectrum is used to compute model function
    bool abs_lifetime_spectrum = true;

public:

    double get_amplitude_threshold(){
        return amplitude_threshold;
    }

    void set_amplitude_threshold(double v){
        amplitude_threshold = v;
    }

    bool get_use_amplitude_threshold(){
        return use_amplitude_threshold;
    }

    void set_use_amplitude_threshold(bool v){
        use_amplitude_threshold = v;
    }

    bool get_abs_lifetime_spectrum() const{
        return abs_lifetime_spectrum;
    }

    void set_abs_lifetime_spectrum(bool v){
        abs_lifetime_spectrum = v;
    }

    void set_lifetime_spectrum(double *input, int n_input) {
        _lifetime_spectrum.resize(n_input);
        _lifetime_spectrum.assign(input, input+n_input);
    }

    void get_lifetime_spectrum(double **output_view, int *n_output) {
        lt = _lifetime_spectrum;
        if(use_amplitude_threshold){
            discriminate_small_amplitudes(
                    lt.data(), lt.size(),
                    amplitude_threshold
            );
        }
        if (abs_lifetime_spectrum) {
            for (int i = 0; i < lt.size(); i++) lt[i] = std::abs(lt[i]);
        }
        *output_view = lt.data();
        *n_output = lt.size();
    }

    DecayLifetimeSpectrum(
            double* lifetime_spectrum = nullptr, int n_lifetime_spectrum = 0,
            bool use_amplitude_threshold = false,
            bool abs_lifetime_spectrum = false,
            double amplitude_threshold = std::numeric_limits<double>::epsilon()
            ){
        set_use_amplitude_threshold(use_amplitude_threshold);
        set_abs_lifetime_spectrum(abs_lifetime_spectrum);
        set_amplitude_threshold(amplitude_threshold);
        set_lifetime_spectrum(lifetime_spectrum, n_lifetime_spectrum);
    }

};


#endif //FIT2X_DECAYLIFETIMESPECTRUM_H
