#ifndef FIT2X_DECAYBACKGROUND_H
#define FIT2X_DECAYBACKGROUND_H

#include <iostream>

# include "DecayModifier.h"

class DecayBackground: public DecayModifier{

private:

    /// The constant offset in the model
    double constant_offset = 0.0;

public:

    void set_constant_offset(double v){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayBackground::set_constant_offset" << std::endl;
#endif

        constant_offset = v;
    }

    double get_constant_offset(){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayBackground::get_constant_offset" << std::endl;
#endif

        return constant_offset;
    }

    void add(DecayCurve& decay) {
        double offset = get_constant_offset();
#ifdef VERBOSE_FIT2X
        std::clog << "DecayBackround::add" << std::endl;
        std::clog << "-- constant_offset:" << offset << std::endl;
#endif
        decay += offset;
    }

    void resize(size_t n){
#ifdef VERBOSE_FIT2X
        std::clog << "DecayBackground::resize" << std::endl;
#endif
        // For future use / non constant background
    }

    DecayBackground(
            double constant_offset = 0.0
            ) {
#ifdef VERBOSE_FIT2X
        std::clog << "DecayBackground::DecayBackground" << std::endl;
#endif

        set_constant_offset(constant_offset);
    }

    ~DecayBackground() = default;

};

#endif //FIT2X_DECAYBACKGROUND_H
