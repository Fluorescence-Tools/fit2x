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
    std::vector<double> _linearization_table = std::vector<double>();

    /// If set to true multiply the linearization to the model function
    int _use_linearization = false;

public:

    void resize(size_t n){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayLinearization::resize" << std::endl;
#endif
        _linearization_table.resize(n, 1.0);
    }

    void set_linearization_table(double *input, int n_input) {
#ifdef VERBOSE_FIT2X
        std::clog << "DecayLinearization::set_linearization_table" << std::endl;
#endif
        if((input != nullptr) && (n_input >= 0)){
            _linearization_table.resize(n_input);
            _linearization_table.assign(input, input + n_input);
        }
    }

    void set_linearization_table(std::vector<double>& v) {
#ifdef VERBOSE_FIT2X
        std::clog << "DecayLinearization::set_linearization_table" << std::endl;
#endif
        _linearization_table = v;
    }

    void get_linearization_table(double **output_view, int *n_output) {
#ifdef VERBOSE_FIT2X
        std::clog << "DecayLinearization::get_linearization_table" << std::endl;
#endif
        *output_view = _linearization_table.data();
        *n_output = _linearization_table.size();
    }

    bool set_use_linearization(bool v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayLinearization::set_use_linearization" << std::endl;
        std::clog << "-- use_linearization:" << v << std::endl;
#endif
        _use_linearization = v;
    }

    bool get_use_linearization(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayLinearization::get_use_linearization" << std::endl;
#endif
        return _use_linearization;
    }

    DecayLinearization(
            double* linearization_table, int n_linearization_table,
            bool use_linearization = false){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayLinearization::DecayLinearization" << std::endl;
#endif
        set_linearization_table(linearization_table, n_linearization_table);
        set_use_linearization(use_linearization);
    }

    void add(DecayCurve* decay){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayLinearization::add" << std::endl;
        std::clog << "-- use_linearization: " << get_use_linearization() << std::endl;
#endif
        if(get_use_linearization()){
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
