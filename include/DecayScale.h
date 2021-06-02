#ifndef FIT2X_DECAYSCALE_H
#define FIT2X_DECAYSCALE_H

#include <iostream>

#include "DecayCurve.h"
#include "DecayModifier.h"
#include "fsconv.h"


class DecayScale : public DecayModifier{

private:

    /// A constant that is subtracted from the data
    double _constant_background = 0.0;

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
            size_t start, size_t stop,
            double constant_background,
            double* model,
            double* data,
            double* squared_data_weights
    ){
#if VERBOSE_FIT2X
        std::clog << "DecayScale::scale_model" << std::endl;
#endif
        double scale = 0.0;
        rescale_w_bg(model, data, squared_data_weights, constant_background, &scale, start, stop);
    }

    /// Number of photons in data between start and stop (if model is scaled to data). Otherwise
    /// user-specified number of photons
    double get_number_of_photons(){
        double re = 0.0;
        DecayCurve* d = get_data();
        for(size_t i = get_start(d); i < get_stop(d); i++)
            re += d->y[i];
        return re;
    }

    double get_constant_background() const {
        return _constant_background;
    }

    void set_constant_background(double v){
        _constant_background = v;
    }

    void set(
            DecayCurve* data = nullptr,
            double constant_background = 0.0,
            int start = 0, int stop = -1, bool active = true
    ){
        set_data(data);
        set_constant_background(constant_background);
        set_start(start);
        set_stop(stop);
        set_active(active);
    }

    DecayScale(
            DecayCurve* data = nullptr,
            double constant_background = 0.0,
            int start=0, int stop=-1, bool active=true
    ) : DecayModifier(data, start, stop, active){
        set_constant_background(constant_background);
    }

    void add(DecayCurve* decay){
        if(is_active()){
            auto d = get_data();
            scale_model(
                get_start(decay),get_stop(decay),
                get_constant_background(),
                decay->y.data(),
                d->y.data(), d->ey.data()
            );
        }
    }

};

#endif //FIT2X_DECAYSCALE_H
