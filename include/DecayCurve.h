#ifndef FIT2X_DECAYCURVE_H
#define FIT2X_DECAYCURVE_H

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream> /* std::cerr */
#include <cstring> /* strcmp */
#include <memory> /* shared_ptr */
#include <limits> /* std::numeric_limits */

#include "tttrlib/TTTR.h"


class DecayCurve: std::enable_shared_from_this<DecayCurve>{

    friend class DecayConvolution;
    friend class DecayPileup;
    friend class DecayLinearization;
    friend class DecayScale;
    friend class DecayScore;
    friend class DecayBackground;
    friend class Decay;

private:

    void compute_weights(const char* noise_model = "poisson");
    int max_validated_size(const DecayCurve& other) const;
    std::vector<double> _y; // original y values
    std::string noise_model = "NA";
    double shift = 0.0;
    bool shift_set_outside = true;
    double shift_outside_value = 0.0;
    /// Acquisition time of the experiment in seconds
    double acquisition_time = 1e9; // to give some room to the top

    void update_weights(){
        w.resize(ey.size());
        w2.resize(ey.size());
        for(int i = 0; i<ey.size();i++){
            if(ey[i] > 0.0){
                w[i] = 1. / (ey[i]);
                w2[i] = 1. / (ey[i] * ey[i]);
            } else{
                w[i] = 0.0;
                w2[i] = 0.0;
            }
        }
    }

    void update_ey(){
        ey.resize(w.size());
        w2.resize(w.size());
        for(int i = 0; i<w.size();i++){
            if(w[i] > 0.0){
                ey[i] = 1. / (w[i]);
                w2[i] = w[i] * w[i];
            } else{
                w2[i] = 0.0;
                ey[i] = std::numeric_limits<double>::infinity();
            }
        }
    }

    std::vector<double> x;
    std::vector<double> y;
    /// error in y
    std::vector<double> ey;
    /// weights, i.e., 1./ey
    std::vector<double> w;
    /// Squared weights, i.e., (1./ey)**2.0
    std::vector<double> w2;

public:

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

    bool restrict_to_positive_values = true;

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
    static void add_arrays(
            double **output_view, int *n_output,
            double *curve1, int n_curve1,
            double *curve2, int n_curve2,
            double areal_fraction_curve2,
            int start = 0,
            int stop = -1
    );

    size_t size() const{
        return x.size();
    }

    bool empty(){
        return x.empty();
    }

    void resize(int n);

    double get_dx();

    std::vector<double> get_x(){
        return x;
    }

    void get_x(double **output_view, int *n_output){
        *output_view = x.data();
        *n_output = x.size();
    }

    void set_x(std::vector<double> v){
        resize(v.size());
        x = v;
    }

    void set_x(double *input, int n_input){
        auto vec = std::vector<double>(input, input+n_input);
        set_x(vec);
    }

    std::vector<double> get_ey(){
        return ey;
    }

    void set_ey(std::vector<double> v){
        ey = v;
        update_weights();
    }

    std::vector<double> get_y(){
        if(restrict_to_positive_values){
            for(auto &d: y){
                d = std::max(0.0, d);
            }
        }
        return y;
    }

    void set_w(std::vector<double> v){
        w = v;
        update_ey();
    }

    void set_w(double* input, int n_input){
        auto vec = std::vector<double>(input, input + n_input);
        set_w(vec);
    }

    std::vector<double> get_w(){
        return w;
    }

    void get_w(double **output_view, int *n_output){
        *output_view = w.data();
        *n_output = w.size();
    }

    void get_y(double **output_view, int *n_output){
        *output_view = y.data();
        *n_output = y.size();
    }


    void set_y(std::vector<double> v){
        resize(v.size());
        y = v; _y = v;
        compute_weights(noise_model.c_str());
    }

    void set_y(double* input, int n_input){
        auto vec = std::vector<double>(input, input + n_input);
        set_y(vec);
    }

    std::vector<double> get_squared_weights(){
        return w2;
    }

    void set_acquisition_time(double v) {
        if(acquisition_time < 0)
            acquisition_time = std::numeric_limits<double>::max();
        acquisition_time = v;
    }

    double get_acquisition_time() const {
        return acquisition_time;
    }

    void set_shift(double v){
        shift = v;
        y = shift_array(_y.data(), _y.size(),
                        shift,
                        shift_set_outside,
                        shift_outside_value);
    }

    bool equal_size(DecayCurve& other){
        return size() == other.size();
    }

    bool equal_size(std::vector<double> other){
        return size() == other.size();
    }

    void set_tttr(
            std::shared_ptr<TTTR> tttr_data,
            int tttr_micro_time_coarsening
    ){
#if VERBOSE_FIT2X
        std::clog << "DecayCurve::set_tttr" << std::endl;
#endif
        if(tttr_data != nullptr){
            double *hist; int n_hist;
            double *time; int n_time;
            TTTR::compute_microtime_histogram(
                    tttr_data.get(),
                    &hist, &n_hist,
                    &time, &n_time,
                    tttr_micro_time_coarsening
            );
            y.assign(hist, hist + n_hist);
            x.assign(time, time + n_time);
            free(hist); free(time);
            auto header = tttr_data->get_header();
            double micro_time_resolution =
                    header->get_micro_time_resolution() / tttr_micro_time_coarsening;
            for(size_t i = 0; i < x.size(); i++) x[i] *= micro_time_resolution;
        }
#if VERBOSE_FIT2X
        else{
            std::clog << "-- no tttr provided" << std::endl;
        }
#endif
    }

    DecayCurve(
            std::vector<double> time_axis,
            std::vector<double> counts,
            double acquisition_time = -1,
            const char* noise_model = "poisson"
            ){
        this->acquisition_time = acquisition_time;
        this->noise_model = noise_model;
        set_x(time_axis);
        set_y(counts);
    }

    DecayCurve(int size = 0, const char* noise_model = "poisson"){
        this->noise_model = noise_model;
        resize(size);
    }

    DecayCurve operator+(const DecayCurve& other) const;
    DecayCurve operator+(double v) const;
    DecayCurve& operator+=(double v);
    DecayCurve operator*(double v) const;
    DecayCurve& operator*=(double v);

};

#endif //FIT2X_DECAYCURVE_H
