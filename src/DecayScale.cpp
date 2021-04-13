#include "DecayScale.h"


void DecayScale::scale_model(
        bool scale_model_to_data,
        double number_of_photons,
        int start, int stop,
        double constant_background,
        double* model,
        double* data,
        double* squared_data_weights
        ){
#if VERBOSE_FIT2X
    std::clog << "-- scaling model to data..." << std::endl;
#endif
    double scale = 1.0;
    if(scale_model_to_data || (number_of_photons < 0)){
#if VERBOSE_FIT2X
        std::clog << "-- constant_background: " << constant_background << std::endl;
        std::clog << "-- scale: " << scale << std::endl;
        std::clog << "-- scaling range: " << start << "," << stop << std::endl;
#endif
        scale = 0.0;
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
