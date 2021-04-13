#ifndef FIT2X_DECAYSCALE_H
#define FIT2X_DECAYSCALE_H

#include "DecayCurve.h"
#include "DecayModifier.h"
#include "fsconv.h"


class DecayScale : public DecayModifier{

private:

    /// A constant that is subtracted from the data
    double _constant_background = 0.0;

    /// Number of photons in output / modified decay.
    double _number_of_photons = -1;

    /// If set to true the model/output is scale to the input data
    bool _scale_model_to_data;

    /// Scaling range (Photons in the decay range (_start, _stop) are used for scaling
    int _scale_start = 0;
    int _scale_stop = -1;

    /// Input data used for scaling
    DecayCurve* data;

public:

    /*!
     * Scales the model in the specified range
     *
     * The model is either scaled to the data or to a specified number of photons
     *
     * @param scale_model_to_data[in] If true the model is scaled to the data
     * @param number_of_photons[in] If scale_model_to_data is false the model is
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

    int get_scale_start(){
        int start = _scale_start;
        int nmax = data->size();
        if(start < 0){
            start = nmax - start;
        }
        return start;
    }

    void set_scale_start(int v){
        _scale_start = v;
    }

    int get_scale_stop(){
        int stop = _scale_stop;
        int nmax = data->size();
        if(stop < 0){
            stop = std::min(nmax, std::max(0, nmax + stop));
        } else if (stop == 0){
            stop = nmax;
        }
        stop = std::min(stop, nmax);
        return stop;
    }

    void set_scale_stop(int v){
        _scale_stop = v;
    }

    /// Number of photons in data between start and stop (if model is scaled to data). Otherwise
    /// user-specified number of photons
    double get_number_of_photons(){
        if(_scale_model_to_data){
            double re = 0.0;
            for(int i=get_scale_start(); i<get_scale_stop(); i++)
                re += data->y[i];
            return re;
        }
        return _number_of_photons;
    }

    void set_number_of_photons(double v){
        _number_of_photons = v;
    }

    void set_scale_model_to_data(bool v){
        _scale_model_to_data = v;
    }

    bool get_scale_model_to_data(){
        return _scale_model_to_data;
    }

    void set_data(DecayCurve* v){
        data = v;
    }

    double get_constant_background(){
        return _constant_background;
    }

    void set_constant_background(double v){
        _constant_background = v;
    }

    void set(
            DecayCurve* data = nullptr,
            bool scale_model_to_data = false,
            double number_of_photons = -1,
            std::vector<int> scale_range = std::vector<int>({0, -1}),
            double constant_background = 0.0
    ){
        set_scale_model_to_data(scale_model_to_data);
        set_data(data);
        set_number_of_photons(number_of_photons);
        set_scale_start(scale_range[0]);
        set_scale_stop(scale_range[1]);
        set_constant_background(constant_background);
    }

    DecayScale(
            DecayCurve* data = nullptr,
            bool scale_model_to_data = false,
            double number_of_photons = -1,
            std::vector<int> scale_range = std::vector<int>({0, -1}),
            double constant_background = 0.0
            ){
        set(
                data,
                scale_model_to_data,
                number_of_photons,
                scale_range,
                constant_background);
    }

    ~DecayScale() = default;

    void resize(size_t n){
        // for potential future use
    }

    void add(DecayCurve* decay){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayScale::add" << std::endl;
#endif
        scale_model(
                get_scale_model_to_data(),
                get_number_of_photons(),
                get_scale_start(),
                get_scale_stop(),
                get_constant_background(),
                decay->y.data(),
                data->y.data(),
                data->get_squared_weights().data()
        );
    }

};

#endif //FIT2X_DECAYSCALE_H
