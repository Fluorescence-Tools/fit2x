#include "DecayCurve.h"


void DecayCurve::add_arrays(
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

std::vector<double> DecayCurve::shift_array(
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
        std::clog << "DecayCurve::shift_array" << std::endl;
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

void DecayCurve::compute_weights(const char* noise_model){
    if(strcmp(noise_model, "poisson") == 0){
        for (int i = 0; i < size(); i++){
            ey[i] = (y[i] <= 0.0) ? 1.0 : std::sqrt(y[i]);
        }
    }
    update_weights();
}

double DecayCurve::get_dx(){
    double dt;
    if(x.size() > 0){
        dt = x[0];
    }
    if(x.size() > 1) {
        dt = x[1] - dt;
    }
    if(dt == 0.0){
        dt = 1.0;
    }
    return dt;
}

void DecayCurve::resize(int n){
    int old_size = size();
#if VERBOSE_FIT2X
    std::clog << "DecayCurve::resize" << std::endl;
    std::clog << "-- old_size:" << old_size << std::endl;
    std::clog << "-- new_size:" << n << std::endl;
#endif

    x.resize(n);
    y.resize(n, 0.0); _y.resize(n, 0.0);
    ey.resize(n, std::numeric_limits<double>::max());
    w.resize(n, 0.0);
    w2.resize(n, 0.0);

    double dx = get_dx();
    for(int i = old_size; i < n; i++){
        if(i > 0) x[i] = x[i - 1] + dx;
    }
}

int DecayCurve::max_validated_size(const DecayCurve& other) const{
    int n = size();
    if(other.size() < n){
        n = other.size();
    }
    return n;
}

DecayCurve DecayCurve::operator+(const DecayCurve& other) const
{
    int n_max = max_validated_size(other);
    auto d = DecayCurve(n_max);
    for(int i=0; i<n_max; i++){
        d.x[i] = x[i];
        d.y[i] = y[i] + other.y[i];
        d.ey[i] = std::sqrt(ey[i]*ey[i] + other.ey[i]*other.ey[i]);
    }
    d.acquisition_time = acquisition_time + other.acquisition_time;
    return d;
}

DecayCurve DecayCurve::operator+(const double v) const
{
    auto d = DecayCurve(size());
    for(int i=0; i<size(); i++){
        d.x[i] = x[i];
        d.y[i] = y[i] + v;
        d.ey[i] = ey[i];
    }
    return d;
}

DecayCurve& DecayCurve::operator+=(const double v)
{
    for(int i=0; i<size(); i++){
        y[i] += v;
    }
    return *this;
}

DecayCurve DecayCurve::operator*(const double v) const
{
    auto d = DecayCurve(size());
    for(int i=0; i<size(); i++){
        d.x[i] = x[i];
        d.y[i] = y[i] * v;
        d.ey[i] = ey[i] * v;
    }
    return d;
}


DecayCurve& DecayCurve::operator*=(const double v)
{
    for(int i=0; i<size(); i++){
        y[i] *= v;
        ey[i] *= v;
    }
    return *this;
}
