#ifndef FIT2X_DECAYLINEARIZATION_H
#define FIT2X_DECAYLINEARIZATION_H

#include <iostream> /* std::cerr */
#include <vector>

#include "DecayCurve.h"
#include "DecayModifier.h"


class DecayLinearization: public DecayModifier{

    friend class DecayCurve;

private:
    /// The model function is multiplied by this vector is _use_linearization is true
    std::vector<double> _linearization_table;

    /// If set to true multiply the linearization to the model function
    bool use_linearization = false;

public:

    void resize(size_t n){
        _linearization_table.resize(n, 1.0);
    }

    void set_linearization_table(double *input, int n_input) {
        _linearization_table.resize(n_input);
        _linearization_table.assign(input, input + n_input);
    }

    void get_linearization_table(double **output_view, int *n_output) {
        *output_view = _linearization_table.data();
        *n_output = _linearization_table.size();
    }

    bool set_use_linearization(bool v){
        use_linearization = v;
    }

    bool get_use_linearization(){
        return use_linearization;
    }

    DecayLinearization(
            std::vector<double> linearization_table  = std::vector<double>(),
            bool use_linearization = false){
        set_linearization_table(linearization_table.data(), linearization_table.size());
        set_use_linearization(use_linearization);
    }

    void add(DecayCurve* decay){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayLinearization::add" << std::endl;
        std::clog << "-- use_linearization: " << use_linearization << std::endl;
#endif
        if(use_linearization){
            if(!decay->equal_size(_linearization_table)){
                std::cerr << "WARNING: Linearization size mismatch." << std::endl;
            }
            int n_max = std::min(_linearization_table.size(), decay->size());
#ifdef VERBOSE_FIT2X
            std::clog << "-- linearization table: ";
#endif
            for(int i = 0; i < n_max; i++){
#ifdef VERBOSE_FIT2X
                std::clog << _linearization_table[i] << " ";
#endif
                decay->y[i] *= _linearization_table[i];
            }
#ifdef VERBOSE_FIT2X
            std::clog << std::endl;
#endif
        }
    }

    ~DecayLinearization() = default;

};

#endif //FIT2X_DECAYLINEARIZATION_H
