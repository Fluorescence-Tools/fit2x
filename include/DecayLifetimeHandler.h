#ifndef FIT2X_DECAYLIFETIMEHANDLER_H
#define FIT2X_DECAYLIFETIMEHANDLER_H

#include <vector>
#include <limits>
#include <iostream>

#include "fsconv.h"

class DecayLifetimeHandler{

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
        return std::abs(amplitude_threshold);
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

    void set_lifetime_spectrum(std::vector<double> v) {
        _lifetime_spectrum = v;
    }

    void add_lifetime(double amplitude, double lifetime) {
        _lifetime_spectrum.emplace_back(amplitude);
        _lifetime_spectrum.emplace_back(lifetime);
    }

    std::vector<double>& get_lifetime_spectrum(){
        lt = _lifetime_spectrum;
        if(use_amplitude_threshold){
            discriminate_small_amplitudes(
                    lt.data(), lt.size(),
                    amplitude_threshold
            );
        }
        if (abs_lifetime_spectrum) {
            for(auto &v: lt) v = std::abs(v);
        }
        return lt;
    }

    void get_lifetime_spectrum(double **output_view, int *n_output) {
        get_lifetime_spectrum();
        *output_view = lt.data();
        *n_output = (int) lt.size();
    }

    DecayLifetimeHandler(
            std::vector<double> lifetime_spectrum = std::vector<double>(),
            bool use_amplitude_threshold = false,
            bool abs_lifetime_spectrum = false,
            double amplitude_threshold = std::numeric_limits<double>::epsilon()
            ){
        set_use_amplitude_threshold(use_amplitude_threshold);
        set_abs_lifetime_spectrum(abs_lifetime_spectrum);
        set_amplitude_threshold(amplitude_threshold);
        set_lifetime_spectrum(lifetime_spectrum);
    }

};


#endif //FIT2X_DECAYLIFETIMEHANDLER_H
